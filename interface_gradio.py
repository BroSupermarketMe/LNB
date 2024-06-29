import gradio as gr
import os
import asyncio
import tkinter as tk
from tkinter.filedialog import askopenfilename, askdirectory, asksaveasfilename, askopenfilenames

def tk_window_asksavefile(init_dir=os.getcwd(), suffix='') -> str:
    window = tk.Tk()
    window.wm_attributes('-topmost', 1)
    window.withdraw()
    filename = asksaveasfilename(initialdir=init_dir, filetypes=[('', suffix)])
    return filename

async def tk_asksavefile_asy(init_dir=os.getcwd(), suffix='') -> str:
    fname = await asyncio.to_thread(tk_window_asksavefile, init_dir, suffix)
    return fname

def my_function(file, lipid_type, lipid_quantity, optional_lipid_type, boundary_condition, save_path):
    # 在这里添加你的处理逻辑
    output_path = os.path.join(save_path, file.name)
    result = f"文件路径：{output_path}, 脂质种类：{lipid_type}, 数量：{lipid_quantity}, 可选脂质种类：{optional_lipid_type}, 边界条件：{boundary_condition}"
    return result

with gr.Blocks() as demo:
    file_input = gr.File(label="上传文件：", file_options=True)
    lipid_type_dropdown = gr.Dropdown(["脂质1", "脂质2"], label="脂质种类：")
    lipid_quantity_slider = gr.Slider(minimum=0, maximum=100, default=50, label="数量：")
    optional_lipid_checkbox = gr.Checkbox("可选脂质种类")
    boundary_condition_dropdown = gr.Dropdown(["条件1", "条件2"], label="边界条件：")
    save_path_button = gr.Button("选择保存路径")
    save_path_button.click(tk_asksavefile_asy, inputs=[], outputs=[file_input])
    submit_button = gr.Button("提交")
    submit_button.click(my_function, inputs=[file_input, lipid_type_dropdown, lipid_quantity_slider,
                                              optional_lipid_checkbox, boundary_condition_dropdown, save_path_button])

    demo.queue(concurrency_count=1022, max_size=2044).launch(server_name="127.0.0.1", inbrowser=True, quiet=True, share=False)
