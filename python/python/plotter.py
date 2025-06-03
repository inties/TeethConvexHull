import matplotlib.pyplot as plt
import numpy as np

class Plotter:
    """
    用于绘制牙齿凸包和样条拟合结果的类。
    """

    def __init__(self, convex_hulls: dict = None, spline_fits: dict = None, overall_hull: np.ndarray = None):
        """
        初始化Plotter。

        Args:
            convex_hulls: 包含凸包结果的字典，格式: {label: np.ndarray}。
            spline_fits: 包含样条拟合结果的字典，格式: {label: np.ndarray}。
        """
        self._convex_hulls = convex_hulls if convex_hulls is not None else {}
        self._spline_fits = spline_fits if spline_fits is not None else {}
        self._overall_hull = overall_hull if overall_hull is not None else []
        self._figure = None
        self._ax = None

    def set_data(self, convex_hulls: dict = None, spline_fits: dict = None):
        """
        设置或更新要绘制的数据。

        Args:
            convex_hulls: 包含凸包结果的字典，格式: {label: np.ndarray}。
            spline_fits: 包含样条拟合结果的字典，格式: {label: np.ndarray}。
        """
        if convex_hulls is not None:
            self._convex_hulls = convex_hulls
        if spline_fits is not None:
            self._spline_fits = spline_fits

    def plot_hulls(self, color='b', linewidth=2):
        """
        绘制凸包。

        Args:
            color: 绘制颜色。
            linewidth: 线条宽度。
        """
        if not self._ax:
             self._figure, self._ax = plt.subplots(figsize=(10, 8))
             self._ax.set_title("Teeth Geometry")
             self._ax.set_xlabel("X")
             self._ax.set_ylabel("Y")
             self._ax.axis('equal')


        print("Plotting convex hulls...")
        for label, hull_pts in self._convex_hulls.items():
            if len(hull_pts) > 0:
                 # 闭合凸包
                 hull_pts_closed = np.vstack([hull_pts, hull_pts[0]])
                 self._ax.plot(hull_pts_closed[:, 0], hull_pts_closed[:, 1],
                              color=color, linewidth=linewidth, label=f'Hull {label}')

    def plot_overall_hull(self, color='b', linewidth=2):
        #hull_pts_closed = np.vstack([self._overall_hull, self._overall_hull[0]])
        if self._overall_hull is not None and not isinstance(self._overall_hull, np.ndarray):
            self._overall_hull = np.array(self._overall_hull)
        self._ax.plot(self._overall_hull[:, 0], self._overall_hull[:, 1],
                    color=color, linewidth=linewidth)
    def plot_splines(self, color='r', linewidth=2):
        """
        绘制样条拟合曲线。

        Args:
            color: 绘制颜色。
            linewidth: 线条宽度。
        """
        if not self._ax:
             self._figure, self._ax = plt.subplots(figsize=(10, 8))
             self._ax.set_title("Teeth Geometry")
             self._ax.set_xlabel("X")
             self._ax.set_ylabel("Y")
             self._ax.axis('equal')

        print("Plotting spline fits...")
        for label, spline_pts in self._spline_fits.items():
             if len(spline_pts) > 0:
                 self._ax.plot(spline_pts[:, 0], spline_pts[:, 1],
                              color=color, linewidth=linewidth, linestyle='--', label=f'Spline {label}')


    def plot_all(self):
        """
        绘制凸包和样条拟合曲线。
        """
        self.plot_hulls(color='b', linewidth=2)
        self.plot_splines(color='r', linewidth=2)
        if self._ax:
             self._ax.legend()


    def show(self):
        """
        显示绘制的图像。
        """
        if self._figure:
             plt.show()
        else:
             print("No plot has been created yet. Call plot_hulls(), plot_splines(), or plot_all() first.")