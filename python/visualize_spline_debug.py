#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
样条曲线调试数据可视化程序
可视化质心、原始样条曲线、更新后的样条曲线和射线
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
import argparse
import os

def load_debug_data(json_path):
    """加载调试JSON数据"""
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"错误：找不到文件 {json_path}")
        return None
    except json.JSONDecodeError as e:
        print(f"错误：JSON解析失败 - {e}")
        return None

def extract_points(point_list):
    """从点列表中提取x和y坐标"""
    if not point_list:
        return [], []
    x_coords = [point['x'] for point in point_list]
    y_coords = [point['y'] for point in point_list]
    return x_coords, y_coords

def plot_spline_debug(data, save_path=None, show_rays=True, show_original=True):
    """绘制样条曲线调试图"""
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    
    # 获取质心
    center = data['center']
    center_x, center_y = center['x'], center['y']
    
    # 绘制质心
    ax.plot(center_x, center_y, 'ro', markersize=10, label='质心', zorder=10)
    ax.add_patch(Circle((center_x, center_y), 0.5, color='red', alpha=0.3, zorder=5))
    
    # 绘制射线（如果启用）
    if show_rays and 'rays' in data:
        rays = data['rays']
        for i, ray in enumerate(rays):
            if i % 5 == 0:  # 只显示每5条射线，避免过于密集
                ax.plot([ray['start_x'], ray['end_x']], 
                       [ray['start_y'], ray['end_y']], 
                       'k--', alpha=0.3, linewidth=0.5, zorder=1)
    
    # 绘制原始样条曲线（如果启用）
    if show_original:
        # 原始内侧样条曲线
        if 'original_inner_spline' in data:
            orig_inner_x, orig_inner_y = extract_points(data['original_inner_spline'])
            if orig_inner_x:
                ax.plot(orig_inner_x, orig_inner_y, 'b-', linewidth=2, alpha=0.7, 
                       label='原始内侧样条曲线', zorder=3)
                ax.scatter(orig_inner_x, orig_inner_y, c='blue', s=20, alpha=0.7, zorder=4)
        
        # 原始外侧样条曲线
        if 'original_outer_spline' in data:
            orig_outer_x, orig_outer_y = extract_points(data['original_outer_spline'])
            if orig_outer_x:
                ax.plot(orig_outer_x, orig_outer_y, 'g-', linewidth=2, alpha=0.7, 
                       label='原始外侧样条曲线', zorder=3)
                ax.scatter(orig_outer_x, orig_outer_y, c='green', s=20, alpha=0.7, zorder=4)
    
    # 绘制更新后的样条曲线
    # 更新后的内侧样条曲线
    if 'updated_inner_spline' in data:
        upd_inner_x, upd_inner_y = extract_points(data['updated_inner_spline'])
        if upd_inner_x:
            ax.plot(upd_inner_x, upd_inner_y, 'b-', linewidth=3, 
                   label='更新后内侧样条曲线', zorder=6)
            ax.scatter(upd_inner_x, upd_inner_y, c='darkblue', s=40, 
                      marker='s', label='更新后内侧点', zorder=7)
    
    # 更新后的外侧样条曲线
    if 'updated_outer_spline' in data:
        upd_outer_x, upd_outer_y = extract_points(data['updated_outer_spline'])
        if upd_outer_x:
            ax.plot(upd_outer_x, upd_outer_y, 'g-', linewidth=3, 
                   label='更新后外侧样条曲线', zorder=6)
            ax.scatter(upd_outer_x, upd_outer_y, c='darkgreen', s=40, 
                      marker='s', label='更新后外侧点', zorder=7)
    
    # 设置图形属性
    ax.set_xlabel('X 坐标', fontsize=12)
    ax.set_ylabel('Y 坐标', fontsize=12)
    ax.set_title('样条曲线角度均匀化调试可视化', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_aspect('equal', adjustable='box')
    
    # 添加统计信息文本
    if 'statistics' in data:
        stats = data['statistics']
        stats_text = f"""统计信息:
原始内侧点数: {stats.get('original_inner_point_count', 'N/A')}
原始外侧点数: {stats.get('original_outer_point_count', 'N/A')}
更新后内侧点数: {stats.get('updated_inner_point_count', 'N/A')}
更新后外侧点数: {stats.get('updated_outer_point_count', 'N/A')}
射线数量: {stats.get('uniform_angle_count', 'N/A')}
最大半径: {stats.get('max_radius', 'N/A'):.2f}"""
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
               verticalalignment='top', bbox=dict(boxstyle='round', 
               facecolor='wheat', alpha=0.8), fontsize=10)
    
    plt.tight_layout()
    
    # 保存图片
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图片已保存到: {save_path}")
    
    plt.show()

def plot_angle_distribution(data, save_path=None):
    """绘制角度分布图"""
    if 'uniform_angles' not in data:
        print("警告：没有找到角度数据")
        return
    
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    angles = data['uniform_angles']
    angles_deg = [np.degrees(angle) for angle in angles]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 角度分布直方图
    ax1.hist(angles_deg, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.set_xlabel('角度 (度)', fontsize=12)
    ax1.set_ylabel('频次', fontsize=12)
    ax1.set_title('角度分布直方图', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # 极坐标图显示角度分布
    ax2 = plt.subplot(122, projection='polar')
    radii = np.ones(len(angles))  # 所有点在同一半径上
    ax2.scatter(angles, radii, c='red', s=30, alpha=0.7)
    ax2.set_ylim(0, 1.2)
    ax2.set_title('角度分布极坐标图', fontsize=14, pad=20)
    
    plt.tight_layout()
    
    if save_path:
        base_name = os.path.splitext(save_path)[0]
        angle_save_path = f"{base_name}_angles.png"
        plt.savefig(angle_save_path, dpi=300, bbox_inches='tight')
        print(f"角度分布图已保存到: {angle_save_path}")
    
    plt.show()

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='可视化样条曲线调试数据')
    parser.add_argument('--json_path', '-j', default='E:/MylabProjects/debug.json',
                       help='调试JSON文件路径')
    parser.add_argument('--save_path', '-s', default=None,
                       help='保存图片的路径')
    parser.add_argument('--no_rays', action='store_true',
                       help='不显示射线')
    parser.add_argument('--no_original', action='store_true',
                       help='不显示原始样条曲线')
    parser.add_argument('--show_angles', action='store_true',
                       help='显示角度分布图')
    
    args = parser.parse_args()
    
    # 加载数据
    data = load_debug_data(args.json_path)
    if data is None:
        return
    
    print("成功加载调试数据")
    print(f"数据包含的键: {list(data.keys())}")
    
    # 绘制主要的样条曲线图
    plot_spline_debug(data, 
                     save_path=args.save_path,
                     show_rays=not args.no_rays,
                     show_original=not args.no_original)
    
    # 如果需要，绘制角度分布图
    if args.show_angles:
        plot_angle_distribution(data, save_path=args.save_path)

if __name__ == "__main__":
    main() 