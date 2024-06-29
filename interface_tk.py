import tkinter as tk
import tkinter.ttk as ttk
import tkinter.filedialog


class AddComponent:
    def __init__(self, type_choices):
        self.frame = tk.Frame()

        self.label_type = tk.Label(self.frame, text="Type:")
        self.label_type.grid(row=0, column=0)

        self.box_type = ttk.Combobox(self.frame)
        self.box_type["value"] = type_choices
        self.box_type.grid(row=0, column=1)

        self.label_ratio = tk.Label(self.frame, text="\tRelative abundance:")
        self.label_ratio.grid(row=0, column=2)

        self.var_ratio = tk.StringVar()
        self.entry_ratio = tk.Entry(self.frame, textvariable=self.var_ratio)
        self.entry_ratio.grid(row=0, column=3)

        self.frame.pack()


class AddLipid(AddComponent):
    lipid_type_choices = ('DPPC', 'DHPC', 'DLPC', 'DMPC', 'DSPC', 'D6PC', 'PSM', 'POPC', 'POP5',
                          'DOPC', 'DAPC', 'DUPC', 'PDPC', 'ODPC', 'DEPC', 'DPPE', 'DHPE', 'DLPE',
                          'DMPE', 'DSPE', 'POPE', 'DOPE', 'PPCS', 'DOPG', 'POPG', 'DOPS', 'POPS',
                          'DPPS', 'DPPG', 'CPG', 'PPG', 'PPT', 'DSMG', 'DSDG', 'DSSQ', 'GMO', 'PI45',
                          'GM1', 'DGDG', 'MGDG', 'SQDG', 'CER', 'GCER', 'DPPI', 'PI', 'PI34',
                          'CDGDG', 'CMGDG', 'CSQDG', 'CSQDB', 'PDGDG', 'PDGDT', 'PMGDG', 'PMGDT',
                          'PSQDG', 'AMA', 'AMA.w', 'KMA', 'MMA', 'CHOL', 'CHOL0', 'CHOL1', 'CHOL2')

    def __init__(self):
        super(AddLipid, self).__init__(AddLipid.lipid_type_choices)


class AddSolvent(AddComponent):
    solvent_type_choices = ('PW', 'IM', 'BMW', 'SPC', 'SPCM', 'FG4W', 'FG4W-MS', 'GLUC', 'FRUC',
                            'SUCR', 'MALT', 'CELL', 'KOJI', 'SOPH', 'NIGE', 'LAMI', 'TREH', 'W',
                            'NA', 'CL', 'Mg', 'K', 'BUT', 'G')

    def __init__(self):
        super(AddSolvent, self).__init__(AddSolvent.solvent_type_choices)


class AddGas(AddComponent):
    gas_type_choices = ('CO2', 'N2', 'O2', 'H2', 'AIR')

    def __init__(self):
        super(AddGas, self).__init__(AddGas.gas_type_choices)


class InterfaceTK:
    def __init__(self):
        # 窗口本身
        self.window = tk.Tk()
        self.window.geometry("800x600")
        self.window.title("LNB generator")

        # ”输入文件列表“文本框
        self.list_input_files = []
        self.text_input_files = tk.Text(self.window, width=100)

        # “添加输入文件”按钮
        self.button_add_input_file = tk.Button(self.window, text="添加输入文件",
                                               command=self.add_input_file)

        # “开始设置参数”按钮
        self.button_start_set_basic = tk.Button(self.window, text="已添加所有输入文件，下一步（设置参数）",
                                                command=self.start_set_basic)

        # “基本参数”布局框
        self.frame_basic = tk.Frame(self.window)

        self.label_d = tk.Label(self.frame_basic,
                                text="Distance between periodic images (nm)")
        self.label_d.grid(row=0, column=0, sticky=tk.E)
        self.var_d = tk.StringVar()
        self.var_d.set('1')
        self.entry_d = tk.Entry(self.frame_basic, textvariable=self.var_d)
        self.entry_d.grid(row=0, column=1)

        self.label_r = tk.Label(self.frame_basic,
                                text="Radius of the Membrane sphere")
        self.label_r.grid(row=1, column=0, sticky=tk.E)
        self.var_r = tk.StringVar()
        self.var_r.set('5')
        self.entry_r = tk.Entry(self.frame_basic, textvariable=self.var_r)
        self.entry_r.grid(row=1, column=1)

        self.label_x = tk.Label(self.frame_basic,
                                text="X dimension or first lattice vector of system (nm)")
        self.label_x.grid(row=2, column=0, sticky=tk.E)
        self.var_x = tk.StringVar()
        self.var_x.set('20')
        self.entry_x = tk.Entry(self.frame_basic, textvariable=self.var_x)
        self.entry_x.grid(row=2, column=1)

        self.label_y = tk.Label(self.frame_basic,
                                text="Y dimension or first lattice vector of system (nm)")
        self.label_y.grid(row=3, column=0, sticky=tk.E)
        self.var_y = tk.StringVar()
        self.var_y.set('20')
        self.entry_y = tk.Entry(self.frame_basic, textvariable=self.var_y)
        self.entry_y.grid(row=3, column=1)

        self.label_z = tk.Label(self.frame_basic,
                                text="Z dimension or first lattice vector of system (nm)")
        self.label_z.grid(row=4, column=0, sticky=tk.E)
        self.var_z = tk.StringVar()
        self.var_z.set('20')
        self.entry_z = tk.Entry(self.frame_basic, textvariable=self.var_z)
        self.entry_z.grid(row=4, column=1)
        
        # ”开始设置脂质“按钮
        self.button_start_set_lipid = tk.Button(self.window, text="下一步（设置脂质）",
                                                command=self.start_set_lipid)

        # “请设置脂质种类“提示
        self.label_set_lipid = tk.Label(self.window,
                                        text="请选择脂质种类，并设置每种脂质的占比（%）。请确保占比之和为100。")

        # 脂质列表
        self.list_lipid = []

        # “添加脂质“按钮
        self.button_add_lipid = tk.Button(self.window, text="添加脂质", command=self.add_lipid)

        # ”设置a“按钮
        self.button_set_a = tk.Button(self.window, text="已添加所有脂质，继续（设置每个脂质的面积）",
                                      command=self.continue_set_a)

        # ”a“布局框
        self.frame_a = tk.Frame(self.window)
        self.label_a = tk.Label(self.frame_a, text="Area per lipid (nm*nm): ")
        self.label_a.grid(row=0, column=0)
        self.var_a = tk.StringVar()
        self.entry_a = tk.Entry(self.frame_a, textvariable=self.var_a)
        self.entry_a.grid(row=0, column=1)

        # “开始设置溶剂”按钮
        self.button_start_set_solvent = tk.Button(self.window, text="下一步（设置溶液）",
                                                  command=self.start_set_solvent)

        # ”请设置溶剂种类“提示
        self.label_set_solvent = tk.Label(self.window,
                                          text="请选择溶剂种类，并设置每种溶剂的占比（%）。请确保占比之和为100。")

        # 溶剂列表
        self.list_solvent = []

        # “添加溶剂“按钮
        self.button_add_solvent = tk.Button(self.window, text="添加溶剂", command=self.add_solvent)

        # ”设置salt“按钮
        self.button_set_salt = tk.Button(self.window, text="已添加完所有溶剂，继续（设置溶液浓度）",
                                         command=self.continue_set_salt)

        # ”salt“布局框
        self.frame_salt = tk.Frame(self.window)
        self.label_salt = tk.Label(self.frame_salt, text="Salt concentration:")
        self.label_salt.grid(row=0, column=0)
        self.var_salt = tk.StringVar()
        self.entry_salt = tk.Entry(self.frame_salt, textvariable=self.var_salt)
        self.entry_salt.grid(row=0, column=1)

        # “开始设置气体”按钮
        self.button_start_set_gas = tk.Button(self.window, text="下一步（设置气体）",
                                              command=self.start_set_gas)

        # “请设置气体种类“提示
        self.label_set_gas = tk.Label(self.window,
                                      text="请选择气体种类，并设置每种气体的占比（%）。请确保占比之和为100。")

        # 气体列表
        self.list_gas = []

        # “添加气体“按钮
        self.button_add_gas = tk.Button(self.window, text="添加气体", command=self.add_gas)

        # ”设置gden“按钮
        self.button_set_gden = tk.Button(self.window, text="已添加完所有气体，继续（设置气体密度）",
                                         command=self.continue_set_gden)

        # ”gden“布局框
        self.frame_gden = tk.Frame(self.window)
        self.label_gden = tk.Label(self.frame_gden, text="Gas density:")
        self.label_gden.grid(row=0, column=0)
        self.var_gden = tk.StringVar()
        self.entry_gden = tk.Entry(self.frame_gden, textvariable=self.var_gden)
        self.entry_gden.grid(row=0, column=1)

        # ”开始设置保存路径“按钮
        self.button_start_set_save_path = tk.Button(self.window, text="下一步（设置保存路径）",
                                                    command=self.start_set_save_path)

        # ”保存路径“文本框
        self.var_save_path = tk.StringVar()
        self.entry_save_path = tk.Entry(self.window, textvariable=self.var_save_path, width=100)

        # “设置保存路径”按钮
        self.button_save_path = tk.Button(self.window, text="设置保存路径", command=self.ask_save_path)
        
        # “确认设置”按钮
        self.button_show_information = tk.Button(self.window, text="下一步（确认设置）",
                                                 command=self.show_all_input)

        # “输出”布局框
        self.frame_show_files = tk.Frame(self.window)
        self.frame_output = tk.Frame(self.window)

        # 输出的属性
        self.output_d = ""
        self.output_r = ""
        self.output_x = ""
        self.output_y = ""
        self.output_z = ""
        self.output_lipid_dict = {}
        self.output_a = ""
        self.output_solvent_dict = {}
        self.output_salt = ""
        self.output_gas_dict = {}
        self.output_gden = ""
        self.output_save_path = ""
        self.output_all = {}

        # ”确认“按钮（退出GUI）
        self.button_finish = tk.Button(self.window, text="确认，开始计算", command=self.finish)

        # 显示窗口，并显示第一部分功能
        self.text_input_files.pack()
        self.button_add_input_file.pack()
        self.button_start_set_basic.pack()

        self.window.mainloop()

        # print("__init__()运行完毕")
        # print(self.output_all)

    def add_input_file(self):
        path = tkinter.filedialog.askopenfilename()
        self.list_input_files.append(path)
        self.text_input_files.insert(tk.END, path+"\n")
        return path

    def start_set_basic(self):
        self.button_add_input_file.pack_forget()
        self.button_start_set_basic.pack_forget()
        self.text_input_files.pack_forget()
        self.frame_basic.pack()
        self.button_start_set_lipid.pack()
        
    def start_set_lipid(self):
        self.frame_basic.pack_forget()
        self.button_start_set_lipid.pack_forget()
        self.label_set_lipid.pack()
        self.button_add_lipid.pack()
        self.button_set_a.pack()

    def add_lipid(self):
        self.button_add_lipid.pack_forget()
        self.button_set_a.pack_forget()
        self.list_lipid.append(AddLipid())
        self.button_add_lipid.pack()
        self.button_set_a.pack()

    def continue_set_a(self):
        self.button_add_lipid.pack_forget()
        self.button_set_a.pack_forget()
        self.frame_a.pack()
        self.button_start_set_solvent.pack()

    def start_set_solvent(self):
        for i in self.list_lipid:
            i.frame.pack_forget()
        self.label_set_lipid.pack_forget()
        self.frame_a.pack_forget()
        self.button_start_set_solvent.pack_forget()
        self.label_set_solvent.pack()
        self.button_add_solvent.pack()
        self.button_set_salt.pack()

    def add_solvent(self):
        self.button_add_solvent.pack_forget()
        self.button_set_salt.pack_forget()
        self.list_solvent.append(AddSolvent())
        self.button_add_solvent.pack()
        self.button_set_salt.pack()

    def continue_set_salt(self):
        self.button_add_solvent.pack_forget()
        self.button_set_salt.pack_forget()
        self.frame_salt.pack()
        self.button_start_set_gas.pack()

    def start_set_gas(self):
        for i in self.list_solvent:
            i.frame.pack_forget()
        self.label_set_solvent.pack_forget()
        self.frame_salt.pack_forget()
        self.button_start_set_gas.pack_forget()
        self.label_set_gas.pack()
        self.button_add_gas.pack()
        self.button_set_gden.pack()

    def add_gas(self):
        self.button_add_gas.pack_forget()
        self.button_set_gden.pack_forget()
        self.list_gas.append(AddGas())
        self.button_add_gas.pack()
        self.button_set_gden.pack()

    def continue_set_gden(self):
        self.button_add_gas.pack_forget()
        self.button_set_gden.pack_forget()
        self.frame_gden.pack()
        self.button_start_set_save_path.pack()

    def start_set_save_path(self):
        for i in self.list_gas:
            i.frame.pack_forget()
        self.label_set_gas.pack_forget()
        self.frame_gden.pack_forget()
        self.button_start_set_save_path.pack_forget()
        self.entry_save_path.pack()
        self.button_save_path.pack()
        self.button_show_information.pack()

    def ask_save_path(self):
        path = tkinter.filedialog.asksaveasfilename()
        self.var_save_path.set(path)
        return path

    def show_all_input(self):
        self.entry_save_path.pack_forget()
        self.button_save_path.pack_forget()
        self.button_show_information.pack_forget()

        self.output_d = self.entry_d.get()
        self.output_r = self.entry_r.get()
        self.output_x = self.entry_x.get()
        self.output_y = self.entry_y.get()
        self.output_z = self.entry_z.get()
        for i in self.list_lipid:
            self.output_lipid_dict[i.box_type.get()] = i.var_ratio.get()
        self.output_a = self.entry_a.get()
        for i in self.list_solvent:
            self.output_solvent_dict[i.box_type.get()] = i.var_ratio.get()
        self.output_salt = self.entry_salt.get()
        for i in self.list_gas:
            self.output_gas_dict[i.box_type.get()] = i.var_ratio.get()
        self.output_gden = self.entry_gden.get()
        self.output_save_path = self.var_save_path.get()

        # 显示信息
        for i in range(len(self.list_input_files)):
            tk.Label(self.frame_show_files, text="输入文件：").grid(row=i, column=0)
            tk.Label(self.frame_show_files, text=self.list_input_files[i])\
                .grid(row=i, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text="Distance between periodic images (nm): ")\
            .grid(row=0, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_d).grid(row=0, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text=f"Radius of the Membrane sphere: ")\
            .grid(row=1, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_r).grid(row=1, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text="X dimension or first lattice vector of system (nm): ")\
            .grid(row=2, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_x).grid(row=2, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text="Y dimension or first lattice vector of system (nm): ")\
            .grid(row=3, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_y).grid(row=3, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text="Z dimension or first lattice vector of system (nm): ")\
            .grid(row=4, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_z).grid(row=4, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text="脂质: ").grid(row=5, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=f"{self.output_lipid_dict}")\
            .grid(row=5, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text="Area per lipid (nm*nm): ")\
            .grid(row=6, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_a).grid(row=6, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text="溶剂: ").grid(row=7, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=f"{self.output_solvent_dict}")\
            .grid(row=7, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text="溶液浓度: ").grid(row=8, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_salt).grid(row=8, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text=f"气体: ").grid(row=9, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=f"{self.output_gas_dict}").grid(row=9, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text=f"气体密度: ").grid(row=10, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_gden).grid(row=10, column=1, sticky=tk.W)

        tk.Label(self.frame_output, text=f"保存路径: ").grid(row=11, column=0, sticky=tk.E)
        tk.Label(self.frame_output, text=self.output_save_path).grid(row=11, column=1, sticky=tk.W)

        self.frame_show_files.pack()
        self.frame_output.pack()
        self.button_finish.pack()

    def finish(self):
        self.output_all['-d'] = self.output_d
        self.output_all['-r'] = self.output_r
        self.output_all['-x'] = self.output_x
        self.output_all['-y'] = self.output_y
        self.output_all['-u'] = self.output_lipid_dict
        self.output_all['-a'] = self.output_a
        self.output_all['-sol'] = self.output_solvent_dict
        self.output_all['-salt'] = self.output_salt
        self.output_all['-gas'] = self.output_gas_dict
        self.output_all['-gden'] = self.output_gden
        self.output_all['-o'] = self.output_save_path
        self.output_all['-f'] = self.list_input_files

        self.window.destroy()


# test_class = InterfaceTK()
