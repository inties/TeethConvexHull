#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OBJ文件顶点过滤器
根据JSON文件中的labels字段过滤OBJ文件中的顶点，移除标签为0的顶点
"""

import json
import os
import sys
from typing import List, Tuple, Dict


class OBJVertexFilter:
    """OBJ文件顶点过滤器类"""
    
    def __init__(self):
        self.vertices = []  # 存储顶点信息 (x, y, z)
        self.faces = []     # 存储面信息
        self.labels = []    # 存储标签信息
        self.other_lines = []  # 存储其他行信息（如材质、法向量等）
    
    def read_obj_file(self, obj_path: str) -> bool:
        """
        读取OBJ文件
        
        Args:
            obj_path: OBJ文件路径
            
        Returns:
            bool: 读取是否成功
        """
        try:
            with open(obj_path, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        self.other_lines.append((line_num, line))
                        continue
                    
                    parts = line.split()
                    if not parts:
                        self.other_lines.append((line_num, line))
                        continue
                    
                    if parts[0] == 'v':  # 顶点
                        if len(parts) >= 4:
                            try:
                                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                                self.vertices.append((x, y, z))
                            except ValueError:
                                print(f"警告：第{line_num}行顶点数据格式错误: {line}")
                                return False
                        else:
                            print(f"警告：第{line_num}行顶点数据不完整: {line}")
                            return False
                    elif parts[0] == 'f':  # 面
                        self.faces.append(line)
                    else:
                        self.other_lines.append((line_num, line))
            
            print(f"成功读取OBJ文件：{len(self.vertices)}个顶点，{len(self.faces)}个面")
            return True
            
        except FileNotFoundError:
            print(f"错误：找不到OBJ文件 {obj_path}")
            return False
        except Exception as e:
            print(f"错误：读取OBJ文件时发生异常 {e}")
            return False
    
    def read_json_labels(self, json_path: str) -> bool:
        """
        读取JSON文件中的labels字段
        
        Args:
            json_path: JSON文件路径
            
        Returns:
            bool: 读取是否成功
        """
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            if 'labels' not in data:
                print("错误：JSON文件中没有找到'labels'字段")
                return False
            
            self.labels = data['labels']
            
            if len(self.labels) != len(self.vertices):
                print(f"警告：标签数量({len(self.labels)})与顶点数量({len(self.vertices)})不匹配")
                print("程序将按照较小的数量进行处理")
            
            print(f"成功读取JSON文件：{len(self.labels)}个标签")
            return True
            
        except FileNotFoundError:
            print(f"错误：找不到JSON文件 {json_path}")
            return False
        except json.JSONDecodeError as e:
            print(f"错误：JSON文件格式错误 {e}")
            return False
        except Exception as e:
            print(f"错误：读取JSON文件时发生异常 {e}")
            return False
    
    def filter_vertices(self) -> Tuple[List[Tuple[float, float, float]], Dict[int, int]]:
        """
        过滤掉标签为0的顶点
        
        Returns:
            Tuple[List, Dict]: 过滤后的顶点列表和顶点索引映射
        """
        filtered_vertices = []
        vertex_mapping = {}  # 原始索引 -> 新索引的映射
        new_index = 0
        
        min_length = min(len(self.vertices), len(self.labels))
        
        for i in range(min_length):
            if self.labels[i] != 0:  # 保留标签不为0的顶点
                filtered_vertices.append(self.vertices[i])
                vertex_mapping[i + 1] = new_index + 1  # OBJ文件索引从1开始
                new_index += 1
        
        print(f"过滤完成：保留了{len(filtered_vertices)}个顶点（移除了{len(self.vertices) - len(filtered_vertices)}个标签为0的顶点）")
        return filtered_vertices, vertex_mapping
    
    def update_faces(self, vertex_mapping: Dict[int, int]) -> List[str]:
        """
        更新面信息，使用新的顶点索引
        
        Args:
            vertex_mapping: 顶点索引映射
            
        Returns:
            List[str]: 更新后的面信息列表
        """
        updated_faces = []
        
        for face_line in self.faces:
            parts = face_line.split()
            if parts[0] != 'f':
                continue
            
            updated_parts = ['f']
            face_valid = True
            
            for vertex_ref in parts[1:]:
                # 处理顶点引用（可能包含纹理坐标和法向量）
                vertex_indices = vertex_ref.split('/')
                try:
                    original_vertex_index = int(vertex_indices[0])
                    if original_vertex_index in vertex_mapping:
                        vertex_indices[0] = str(vertex_mapping[original_vertex_index])
                        updated_parts.append('/'.join(vertex_indices))
                    else:
                        # 该顶点已被过滤掉，跳过这个面
                        face_valid = False
                        break
                except ValueError:
                    print(f"警告：面定义中的顶点索引格式错误: {vertex_ref}")
                    face_valid = False
                    break
            
            if face_valid and len(updated_parts) >= 4:  # 至少需要3个顶点构成面
                updated_faces.append(' '.join(updated_parts))
        
        print(f"面信息更新完成：保留了{len(updated_faces)}个面（原来有{len(self.faces)}个面）")
        return updated_faces
    
    def write_obj_file(self, output_path: str, filtered_vertices: List[Tuple[float, float, float]], 
                      updated_faces: List[str]) -> bool:
        """
        写入新的OBJ文件
        
        Args:
            output_path: 输出文件路径
            filtered_vertices: 过滤后的顶点列表
            updated_faces: 更新后的面列表
            
        Returns:
            bool: 写入是否成功
        """
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                # 写入文件头注释
                f.write("# OBJ文件 - 已过滤标签为0的顶点\n")
                f.write(f"# 顶点数量: {len(filtered_vertices)}\n")
                f.write(f"# 面数量: {len(updated_faces)}\n\n")
                
                # 写入顶点
                for vertex in filtered_vertices:
                    f.write(f"v {vertex[0]:.6f} {vertex[1]:.6f} {vertex[2]:.6f}\n")
                
                f.write("\n")
                
                # 写入面
                for face in updated_faces:
                    f.write(f"{face}\n")
                
                # 写入其他信息（如果需要的话）
                if self.other_lines:
                    f.write("\n# 其他信息\n")
                    for _, line in self.other_lines:
                        if not line.startswith('v ') and not line.startswith('f '):
                            f.write(f"{line}\n")
            
            print(f"成功保存过滤后的OBJ文件: {output_path}")
            return True
            
        except Exception as e:
            print(f"错误：写入OBJ文件时发生异常 {e}")
            return False
    
    def process(self, obj_path: str, json_path: str, output_path: str) -> bool:
        """
        执行完整的处理流程
        
        Args:
            obj_path: 输入OBJ文件路径
            json_path: 输入JSON文件路径
            output_path: 输出OBJ文件路径
            
        Returns:
            bool: 处理是否成功
        """
        print("开始处理OBJ文件和JSON文件...")
        
        # 读取OBJ文件
        if not self.read_obj_file(obj_path):
            return False
        
        # 读取JSON文件
        if not self.read_json_labels(json_path):
            return False
        
        # 过滤顶点
        filtered_vertices, vertex_mapping = self.filter_vertices()
        
        if not filtered_vertices:
            print("警告：没有保留任何顶点！")
            return False
        
        # 更新面信息
        updated_faces = self.update_faces(vertex_mapping)
        
        # 写入新的OBJ文件
        return self.write_obj_file(output_path, filtered_vertices, updated_faces)


def main():
    """主函数"""
    if len(sys.argv) != 4:
        print("用法: python filter_obj_vertices.py <输入OBJ文件> <输入JSON文件> <输出OBJ文件>")
        print("示例: python filter_obj_vertices.py input.obj labels.json output.obj")
        sys.exit(1)
    
    obj_path = sys.argv[1]
    json_path = sys.argv[2]
    output_path = sys.argv[3]
    
    # 检查输入文件是否存在
    if not os.path.exists(obj_path):
        print(f"错误：OBJ文件不存在: {obj_path}")
        sys.exit(1)
    
    if not os.path.exists(json_path):
        print(f"错误：JSON文件不存在: {json_path}")
        sys.exit(1)
    
    # 创建输出目录（如果不存在）
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 执行处理
    filter_obj = OBJVertexFilter()
    if filter_obj.process(obj_path, json_path, output_path):
        print("\n处理完成！")
    else:
        print("\n处理失败！")
        sys.exit(1)


if __name__ == "__main__":
    main() 