from  information.usrdefined import user_define
from interface_cli import InterfaceCLI
from protein import protein
from membrane import membrane
from solvent import solvent
from output import output

version = "20211229.LNB"

# This script was written based on the open source script insane.py (Djurre H. de Jong et al. J. Chem. Theory Comput. 2013, 9, 687-697.) by Yuan He and Yuxuan Wang in Dr. Xubo Lin's group@Beihang University (https://linxubo.github.io).
# This script can be used to generate lipid nanobubble (spherical lipid monolayer) using one-line command such as "python3 LNB.py -d 1 -r 5 -x 20 -y 20 -z 20 -u DPPC:50 -u DPPE:50 -sol W -salt 0.15 -a 1 -gas N2 -gden 200 -o output.gro".
# 🔍示例参数："-d 1 -r 5 -x 20 -y 20 -z 20 -u DPPC:50 -u DPPE:50 -sol W -salt 0.15 -a 1 -gas N2 -gden 200 -o output_debug.gro"
# 1️⃣目前界面需要设置参数：-d, -r, -x, -y, -z, -u, -sol, -salt, -a, -gas, -gden, -o
# 2️⃣-d -r -x -y -z 给定默认值（数值按示例设置）
# 3️⃣-u -sol -gas 下拉框（具体样式可以参照interface_gradio，-u需要有数值条，具体选项在information/predefined.py）
#   下拉框里面选项的详情在information/predefined.py里面
#   在终端用命令“python .\interface_gradio.py”，显示gradio界面，可以看到具体的下拉框
# 4️⃣-a -salt -gden 分别在用户手动设定-u -salt -gden后设置
# 5️⃣-o默认./output.gro (linux/MacOS)或.\output.gro (Windows) 用户可以手动设置
# Details about the parameters for the one-line command can be found using "python3 main.py -h".
# ⬆️如果你在系统中通过微软应用商城预先安装了非conda通道的Python，请使用"python main.py -h"。
# This script can only be used with python 3

# 交互界面给core的输出：interface_args = ["-d", 1, "-r", 5 -x 20 -y 20 -z 20 -u DPPC:50 -u DPPE:50 -sol W -salt 0.15 -a 1 -gas N2 -gden 200 -o output_debug.gro]


if __name__ == '__main__':
    interface = InterfaceCLI()
    variables, options = interface.get_args()
    arg_list = user_define(variables, options)
    arg_list = protein(*arg_list)
    arg_list = membrane(*arg_list)
    arg_list = solvent(*arg_list)
    output(*arg_list)