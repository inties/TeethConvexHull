###用于可视化凸包求轴嵴连线结果的程序
import json
import matplotlib.pyplot as plt
import numpy as np

def load_json(filepath):
    with open(filepath, 'r', encoding='utf-8') as file:
        return json.load(file)

def plot_points(data,idx):
    plt.figure(figsize=(10, 8))

    for key,points in data.get("resultpoints").items():
        x_vals = [p["x"] for p in points]
        y_vals = [p["y"] for p in points]
        plt.scatter(x_vals, y_vals, label=f'points {key}', s=1)
    if "debugPoints" in data:
        x_vals = [p["x"] for p in data["debugPoints"]]
        y_vals = [p["y"] for p in data["debugPoints"]]
        plt.plot(x_vals, y_vals, 'g--', label='Inner Spline')
    # # 绘制 pointsMap
    # for key, points in data.get(f"pointsMap_{idx}", {}).items():
    #     x_vals = [p["x"] for p in points]
    #     y_vals = [p["y"] for p in points]
    #     plt.scatter(x_vals, y_vals, label=f'pointsMap {key}', s=1)
    
    # # # 绘制 inner spline
    # if f"spline_out_{idx}" in data:
    #     x_vals = [p["x"] for p in data[f"spline_out_{idx}"]]
    #     y_vals = [p["y"] for p in data[f"spline_out_{idx}"]]
    #     plt.plot(x_vals, y_vals, 'g--', label='Inner Spline')
    
    # # 绘制 outer spline
    # if f"spline_out_{idx}" in data:
    #     x_vals = [p["x"] for p in data[f"spline_out_{idx}"]]
    #     y_vals = [p["y"] for p in data[f"spline_out_{idx}"]]
    #     plt.plot(x_vals, y_vals, 'm--', label='Outer Spline')
    # if f"hull_in_{idx}" in data:
    #     x_vals = [p["x"] for p in data[f"hull_in_{idx}"]]
    #     y_vals = [p["y"] for p in data[f"hull_in_{idx}"]]
    #     plt.scatter(x_vals, y_vals, label='Hull In',s=5,c='r')
    # if f"hull_out_{idx}" in data:
    #     x_vals = [p["x"] for p in data[f"hull_out_{idx}"]]
    #     y_vals = [p["y"] for p in data[f"hull_out_{idx}"]]
    #     plt.scatter(x_vals, y_vals, label='Hull Out',s=5,c='b')
    # # 绘制带标签的innerSplineMap
    # if "innerSplineMap" in data:
    #     # 创建颜色映射
    #     colors = plt.cm.rainbow(np.linspace(0, 1, len(data["innerSplineMap"])))
    #     color_dict = {key: colors[i % len(colors)] for i, key in enumerate(data["innerSplineMap"].keys())}
        
    #     for key, points in data["innerSplineMap"].items():
    #         x_vals = [p["x"] for p in points]
    #         y_vals = [p["y"] for p in points]
    #         plt.scatter(x_vals, y_vals, color=color_dict[key], label=f'Inner Spline {key}', s=30, marker='o')
    
    # # 绘制带标签的outerSplineMap
    # if "outerSplineMap" in data:
    #     # 创建颜色映射
    #     colors = plt.cm.rainbow(np.linspace(0, 1, len(data["outerSplineMap"])))
    #     color_dict = {key: colors[i % len(colors)] for i, key in enumerate(data["outerSplineMap"].keys())}
        
    #     for key, points in data["outerSplineMap"].items():
    #         x_vals = [p["x"] for p in points]
    #         y_vals = [p["y"] for p in points]
    #         plt.scatter(x_vals, y_vals, color=color_dict[key], label=f'Outer Spline {key}', s=30, marker='s')
    
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Visualization of JSON Data")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    filepath = r"E:\MylabProjects\debugPoints.json"# 替换为你的 JSON 文件路径
    data = load_json(filepath)
    plot_points(data,0)
    
 
    
