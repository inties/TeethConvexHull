#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
简化版样条曲线可视化程序
只绘制spline_in_X和spline_out_X数据
"""

import json
import matplotlib.pyplot as plt
import numpy as np
import argparse
import colorsys

def load_data(json_path):
    """加载JSON数据"""
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"错误：无法加载文件 {json_path} - {e}")
        return None

def extract_spline_data(data):
    """提取样条曲线数据"""
    splines = {}
    
    # 提取所有spline_in_X和spline_out_X数据
    for key, value in data.items():
        if key.startswith('spline_in_') and key.split('_')[-1].isdigit():
            idx = int(key.split('_')[-1])
            if idx not in splines:
                splines[idx] = {}
            splines[idx]['inner'] = value
        elif key.startswith('spline_out_') and key.split('_')[-1].isdigit():
            idx = int(key.split('_')[-1])
            if idx not in splines:
                splines[idx] = {}
            splines[idx]['outer'] = value
    
    return splines

def extract_points(point_list):
    """从点列表提取x,y坐标"""
    if not point_list:
        return [], []
    return [p['x'] for p in point_list], [p['y'] for p in point_list]

def generate_colors(n):
    """生成n种不同颜色"""
    return [colorsys.hsv_to_rgb(i/n, 0.8, 0.9) for i in range(n)]

def plot_splines(splines, save_path=None):
    """绘制样条曲线"""
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    fig, ax = plt.subplots(figsize=(12, 10))
    colors = generate_colors(len(splines))
    
    for i, (idx, spline_data) in enumerate(sorted(splines.items())):
        if idx != 3:  # 只绘制索引为4的曲线
            continue
        
        color = colors[i]
        # 绘制内侧样条曲线
        if 'inner' in spline_data:
            x, y = extract_points(spline_data['inner'])
            if x:
                ax.plot(x, y, '-', color=color, linewidth=2.5, 
                       label=f'内侧样条 {idx}', alpha=0.8)
                ax.scatter(x, y, c=[color], s=25, marker='o', alpha=0.7)
        
        # 绘制外侧样条曲线
        if 'outer' in spline_data:
            x, y = extract_points(spline_data['outer'])
            if x:
                ax.plot(x, y, '--', color=color, linewidth=2.5, 
                       label=f'外侧样条 {idx}', alpha=0.8)
                ax.scatter(x, y, c=[color], s=25, marker='s', alpha=0.7)
    
    # 设置图形属性
    ax.set_xlabel('X 坐标', fontsize=12)
    ax.set_ylabel('Y 坐标', fontsize=12)
    ax.set_title(f'样条曲线可视化 (共{len(splines)}组)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存到: {save_path}")
    
    plt.show()

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='简化版样条曲线可视化')
    parser.add_argument('--json_path', '-j', default='E:/MylabProjects/debug.json',
                       help='JSON文件路径')
    parser.add_argument('--save_path', '-s', default=None,
                       help='保存图片路径')
    
    args = parser.parse_args()
    
    # 加载和处理数据
    data = load_data(args.json_path)
    if data is None:
        return
    
    print(f"加载数据成功，包含键: {len(data)} 个")
    
    splines = extract_spline_data(data)
    if not splines:
        print("错误：未找到spline_in_X或spline_out_X数据")
        return
    
    print(f"找到 {len(splines)} 组样条曲线数据")
    
    # 绘制图形
    plot_splines(splines, args.save_path)

if __name__ == "__main__":
    main()