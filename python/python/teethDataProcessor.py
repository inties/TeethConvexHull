import numpy as np
# 导入cv2和可能的样条库，但在接口中不实现
import cv2
# from scipy.interpolate import splrep, splev # Example spline library

class DataProcessor:
    """
    用于处理牙齿点数据的类，计算凸包和样条拟合。
    """

    def __init__(self, points_map: dict):
        """
        初始化DataProcessor。

        Args:
            points_map: 包含牙齿点数据的字典。
                        例如：{ "1": [{"x": 100, "y": 200}, ...], "2": [...], ... }
        """
        self._points_map = points_map
        self._convex_hulls = {}  # 存储凸包结果，格式: {label: np.ndarray}
        self._spline_fits = {}   # 存储样条拟合结果，格式: {label: np.ndarray}
        self.calculate_convex_hulls()
        self.calculate_overall_hull()

    def calculate_convex_hulls(self):
        """
        计算points_map中每个牙齿的凸包。
        结果存储在内部的 _convex_hulls 属性中。
        """
        # TODO: 实现凸包计算逻辑，使用 cv2.convexHull 或其他方法
        print("Calculating convex hulls...")
        # Example placeholder:
        for label, points in self._points_map.items():
            pts = np.array([[p["x"], p["y"]] for p in points], dtype=np.float32)
            if len(pts) >= 3:
                hull_indices = cv2.convexHull(pts, returnPoints=False)
                self._convex_hulls[label] = pts[hull_indices[:, 0]]
    def calculate_overall_hull(self):
        """
        使用所有已计算的单个牙齿凸包点，计算整个牙齿区域的整体凸包。
        结果存储在内部的 _overall_hull 属性中。
        """
        if not self._convex_hulls:
            print("Cannot calculate overall hull: Individual convex hulls have not been calculated yet or are empty.")
            self._overall_hull = None
            return

        print("Calculating overall convex hull from individual hull points...")

        # 收集所有单个牙齿凸包的顶点到一个列表中
        all_hull_points_list = []
        for label, hull_pts in self._convex_hulls.items():
            if hull_pts.shape[0] > 0: # 只添加非空的凸包点
                 all_hull_points_list.append(hull_pts)

        # 检查是否有足够的点来计算整体凸包
        if not all_hull_points_list:
             print("No points available from individual hulls to calculate overall hull.")
             self._overall_hull = None
             return

        # 将所有点列表合并成一个大的numpy数组
        all_hull_points = np.vstack(all_hull_points_list)

        # 计算整体凸包
        # 确保总点数 >= 3
        if len(all_hull_points) >= 3:

            # 使用 returnPoints=True 直接获取凸包顶点坐标
            self._overall_hull = cv2.convexHull(all_hull_points, returnPoints=True)
            # cv2.convexHull(returnPoints=True) 返回的形状是 (Q, 1, 2)。
            # 我们通常需要形状 (Q, 2)，所以这里进行 reshape
            hull_indices = cv2.convexHull(all_hull_points, returnPoints=False)
            self._overall_hull = all_hull_points[hull_indices[:, 0]]
            ##列表转numpy
            self._overall_hull = np.array(self._overall_hull)

        #     except cv2.error as e:
        #          print(f"Error calculating overall convex hull: {e}")
        #          self._overall_hull = np.array([], dtype=np.float32).reshape(-1, 2)

        # else:
        #     print(f"Not enough points ({len(all_hull_points)}) from individual hulls to calculate overall hull.")
        #     self._overall_hull = np.array([], dtype=np.float32).reshape(-1, 2) # Store empty array if not possible
    def calculate_spline_fits(self, num_points: int = 100):
        """
        对points_map中每个牙齿的点进行样条拟合。

        Args:
            num_points: 样条拟合后生成的点数量。
        结果存储在内部的 _spline_fits 属性中。
        """
        # TODO: 实现样条拟合逻辑，例如使用 scipy.interpolate
        print(f"Calculating spline fits with {num_points} points...")
        # Example placeholder:
        # for label, points in self._points_map.items():
        #     pts = np.array([[p["x"], p["y"]] for p in points], dtype=np.float64)
        #     # Ensure points are ordered for spline fitting (e.g., by angle from centroid)
        #     # ... ordering logic ...
        #     if len(pts) >= 4: # Spline usually needs more points than hull
        #         tck, u = splrep(pts[:, 0], pts[:, 1], k=3) # Example B-spline of degree 3
        #         u_new = np.linspace(0, 1, num_points)
        #         x_new, y_new = splev(u_new, tck)
        #         self._spline_fits[label] = np.vstack([x_new, y_new]).T
        pass # Implementation goes here
    def get_convex_hulls(self) -> dict:
        return self._convex_hulls
    def get_overall_hull(self):
        if self._overall_hull is not None and not isinstance(self._overall_hull, np.ndarray):
            self._overall_hull = np.array(self._overall_hull)
        return self._overall_hull

    def get_spline_fits(self) -> dict:
        """
        获取计算得到的样条拟合结果。

        Returns:
            一个字典，键是牙齿标签 (str)，值是表示样条曲线点的numpy数组 (np.ndarray)。
            如果尚未计算，则返回空字典。
        """
        return self._spline_fits