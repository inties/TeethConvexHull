#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OBJ文件顶点过滤器 - 交互式版本
根据JSON文件中的labels字段过滤OBJ文件中的顶点，移除标签为0的顶点
"""

import json
import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from filter_obj_vertices import OBJVertexFilter


class OBJFilterGUI:
    """OBJ文件过滤器图形界面"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("OBJ文件顶点过滤器")
        self.root.geometry("600x400")
        
        self.obj_path = tk.StringVar()
        self.json_path = tk.StringVar()
        self.output_path = tk.StringVar()
        
        self.setup_ui()
    
    def setup_ui(self):
        """设置用户界面"""
        # 主框架
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 标题
        title_label = ttk.Label(main_frame, text="OBJ文件顶点过滤器", font=("Arial", 16, "bold"))
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # OBJ文件选择
        ttk.Label(main_frame, text="OBJ文件:").grid(row=1, column=0, sticky=tk.W, pady=5)
        ttk.Entry(main_frame, textvariable=self.obj_path, width=50).grid(row=1, column=1, padx=5, pady=5)
        ttk.Button(main_frame, text="浏览", command=self.browse_obj_file).grid(row=1, column=2, pady=5)
        
        # JSON文件选择
        ttk.Label(main_frame, text="JSON文件:").grid(row=2, column=0, sticky=tk.W, pady=5)
        ttk.Entry(main_frame, textvariable=self.json_path, width=50).grid(row=2, column=1, padx=5, pady=5)
        ttk.Button(main_frame, text="浏览", command=self.browse_json_file).grid(row=2, column=2, pady=5)
        
        # 输出文件设置
        ttk.Label(main_frame, text="输出文件:").grid(row=3, column=0, sticky=tk.W, pady=5)
        ttk.Entry(main_frame, textvariable=self.output_path, width=50).grid(row=3, column=1, padx=5, pady=5)
        ttk.Button(main_frame, text="保存为", command=self.browse_output_file).grid(row=3, column=2, pady=5)
        
        # 处理按钮
        process_button = ttk.Button(main_frame, text="开始处理", command=self.process_files)
        process_button.grid(row=4, column=1, pady=20)
        
        # 进度条
        self.progress = ttk.Progressbar(main_frame, mode='indeterminate')
        self.progress.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        # 日志显示区域
        ttk.Label(main_frame, text="处理日志:").grid(row=6, column=0, sticky=tk.W, pady=(20, 5))
        
        log_frame = ttk.Frame(main_frame)
        log_frame.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        self.log_text = tk.Text(log_frame, height=12, wrap=tk.WORD)
        scrollbar = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        # 配置网格权重
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(7, weight=1)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
    
    def browse_obj_file(self):
        """选择OBJ文件"""
        filename = filedialog.askopenfilename(
            title="选择OBJ文件",
            filetypes=[("OBJ files", "*.obj"), ("All files", "*.*")]
        )
        if filename:
            self.obj_path.set(filename)
            # 自动设置输出文件名
            if not self.output_path.get():
                base_name = os.path.splitext(filename)[0]
                self.output_path.set(f"{base_name}_filtered.obj")
    
    def browse_json_file(self):
        """选择JSON文件"""
        filename = filedialog.askopenfilename(
            title="选择JSON文件",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if filename:
            self.json_path.set(filename)
    
    def browse_output_file(self):
        """选择输出文件"""
        filename = filedialog.asksaveasfilename(
            title="保存过滤后的OBJ文件",
            filetypes=[("OBJ files", "*.obj"), ("All files", "*.*")],
            defaultextension=".obj"
        )
        if filename:
            self.output_path.set(filename)
    
    def log_message(self, message):
        """添加日志消息"""
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.root.update()
    
    def process_files(self):
        """处理文件"""
        # 验证输入
        if not self.obj_path.get():
            messagebox.showerror("错误", "请选择OBJ文件")
            return
        
        if not self.json_path.get():
            messagebox.showerror("错误", "请选择JSON文件")
            return
        
        if not self.output_path.get():
            messagebox.showerror("错误", "请设置输出文件路径")
            return
        
        # 清空日志
        self.log_text.delete(1.0, tk.END)
        
        # 开始进度条
        self.progress.start()
        
        try:
            # 创建过滤器并重定向输出到日志
            filter_obj = OBJVertexFilter()
            
            # 重写print函数以便输出到GUI
            original_print = print
            def gui_print(*args, **kwargs):
                message = ' '.join(str(arg) for arg in args)
                self.log_message(message)
                original_print(*args, **kwargs)
            
            # 临时替换print函数
            import builtins
            builtins.print = gui_print
            
            # 执行处理
            success = filter_obj.process(
                self.obj_path.get(),
                self.json_path.get(),
                self.output_path.get()
            )
            
            # 恢复print函数
            builtins.print = original_print
            
            if success:
                self.log_message("\n✅ 处理完成！")
                messagebox.showinfo("成功", "文件处理完成！")
            else:
                self.log_message("\n❌ 处理失败！")
                messagebox.showerror("失败", "文件处理失败，请查看日志了解详情。")
        
        except Exception as e:
            self.log_message(f"\n❌ 发生异常: {e}")
            messagebox.showerror("错误", f"处理过程中发生异常：{e}")
        
        finally:
            # 停止进度条
            self.progress.stop()


def main():
    """主函数"""
    root = tk.Tk()
    app = OBJFilterGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main() 