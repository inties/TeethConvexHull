#计算每颗牙齿的凸包
import cv2
import numpy as np
import json
import matplotlib.pyplot as plt
from teethDataProcessor import DataProcessor
from plotter import Plotter
##pointsmap
##

# 1. 读取json文件
with open(r"E:\MylabProjects\ConvexHull\convex\data\OEPN5LUR\OEPN5LUR_lower_out.json", 'r', encoding='utf-8') as f:
    data = json.load(f)

# 假设json结构为 { "pointsMap": { "1": [{"x":..., "y":...}, ...], "2": [...], ... } }
points_map = data["pointsMap_0"]
teeth_data_processor = DataProcessor(points_map)
teeth_plotter = Plotter(teeth_data_processor.get_convex_hulls(),teeth_data_processor.get_overall_hull)
teeth_plotter.plot_hulls()
teeth_plotter.plot_overall_hull()
teeth_plotter.show()

