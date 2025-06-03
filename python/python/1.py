import open3d as o3d
import numpy as np
import sys
import os

def load_and_visualize_lines(obj_path):
    """
    加载 .obj 文件中的线段并使用 Open3D 可视化
    """
    if not os.path.exists(obj_path):
        print(f"错误: 文件 {obj_path} 不存在")
        return
    
    # 手动读取 .obj 文件
    vertices = []
    lines = []
    with open(obj_path, 'r') as f:
        for line in f:
            if line.startswith('v '):
                parts = line.split()
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif line.startswith('l '):
                parts = line.split()
                # 转换为 0-based 索引
                line_indices = [int(p) - 1 for p in parts[1:]]
                lines.append(line_indices)
    
    # 创建 LineSet 对象
    line_set = o3d.geometry.LineSet()
    line_set.points = o3d.utility.Vector3dVector(np.array(vertices))
    line_set.lines = o3d.utility.Vector2iVector(np.array(lines))
    # 设置线条颜色为红色 (RGB值范围为0-1)
    line_set.colors = o3d.utility.Vector3dVector(np.tile([1.0, 0.0, 0.0], (len(lines), 1)))
    
    # 打印信息
    print(f"线框包含 {len(line_set.points)} 个顶点和 {len(line_set.lines)} 条线段")
    
    # 创建可视化窗口
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    
    # 添加线框到可视化器
    vis.add_geometry(line_set)
    
    # 设置渲染选项
    opt = vis.get_render_option()
    opt.background_color = np.asarray([0.1, 0.1, 0.1])  # 深灰色背景
    
    # 运行可视化
    vis.run()
    vis.destroy_window()

def load_and_visualize_lines_utf8(obj_path):
    if not os.path.exists(obj_path):
        print(f"错误: 文件 {obj_path} 不存在")
        return
    
    # 手动读取 .obj 文件，指定 UTF-8 编码
    vertices = []
    lines = []
    try:
        with open(obj_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('v '):
                    parts = line.split()
                    vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
                elif line.startswith('l '):
                    parts = line.split()
                    # 转换为 0-based 索引
                    line_indices = [int(p) - 1 for p in parts[1:]]
                    lines.append(line_indices)
    except UnicodeDecodeError as e:
        print(f"编码错误: {e}")
        print("尝试使用其他编码（例如 'latin1' 或 'gbk'）打开文件")
        return
    
    # 创建 LineSet 对象
    line_set = o3d.geometry.LineSet()
    line_set.points = o3d.utility.Vector3dVector(np.array(vertices))
    line_set.lines = o3d.utility.Vector2iVector(np.array(lines))
    # 设置线条颜色为红色 (RGB值范围为0-1)
    line_set.colors = o3d.utility.Vector3dVector(np.tile([1.0, 0.0, 0.0], (len(lines), 1)))
    
    # 打印信息
    print(f"线框包含 {len(line_set.points)} 个顶点和 {len(line_set.lines)} 条线段")
    
    # 创建可视化窗口
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    
    # 添加线框到可视化器
    vis.add_geometry(line_set)
    
    # 设置渲染选项
    opt = vis.get_render_option()
    opt.background_color = np.asarray([0.1, 0.1, 0.1])  # 深灰色背景
    
    # 运行可视化
    vis.run()
    vis.destroy_window()

if __name__ == "__main__":
    # 默认的 obj 文件路径
    default_obj_path = r"E:\\MylabProjects\\ConvexHull\\convex\\convex1\\testModels\\teethTest_lower.obj"
    if len(sys.argv) == 2:
        obj_path = sys.argv[1]
    else:
        print(f"未提供 obj 文件路径，使用默认路径: {default_obj_path}")
        obj_path = default_obj_path
    load_and_visualize_lines(obj_path)