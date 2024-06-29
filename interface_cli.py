import sys
from utils.structure import Structure
from utils.algebra import vector
from utils.parse_command import Option
from interface_tk import InterfaceTK

class InterfaceCLI():
    def __init__(self):
        self.tm = []
        self.lipL = []
        self.lipU = []
        self.solv = []
        self.gas = []

        # HII edit - lipid definition, for extra lipid definitaions
        self.usrmols = []
        self.usrheads = []
        self.usrlinks = []
        self.usrtails = []
        self.usrLipHeadMapp = {  # Define supported lipid head beads. One letter name mapped to atom name
            "C": ('NC3'),  # NC3 = Choline
            "E": ('NH3'),  # NH3 = Ethanolamine
            "G": ('GL0'),  # GL0 = Glycerol
            "S": ('CNO'),  # CNO = Serine
            "P": ('PO4')  # PO4 = Phosphate
        }
        self.usrIndexToLetter = "A B C D E F G H I J K L M N".split()  # For naming lipid tail beads

        # Description
        self.desc = "The following are the functions that the program can perform:"
        # Option list
        self.options = [
            #   option           type number default description
            # HII edit - lipid definition (last options are for additional lipid specification)
            ("""
        Input/output related options
        """),
            ("-f", Option(self.tm.append, 1, None, "Input GRO or PDB file 1: Protein")),
            ("-o", Option(str, 1, None, "Output GRO file: Membrane with Protein")),
            ("-p", Option(str, 1, None, "Optional rudimentary topology file")),
            """
        Periodic boundary conditions 
        If -d is given, set up PBC according to -pbc such that no periodic
        images are closer than the value given.  This will make the numbers
        provided for lipids be interpreted as relative numbers. If -d is
        omitted, those numbers are interpreted as absolute numbers, and the
        PBC are set to fit the given number of lipids in.
        """,
            ("-pbc", Option(str, 1, "hexagonal", "PBC type: hexagonal, rectangular, square, cubic, optimal or keep")),
            ("-d", Option(float, 1, None, "Distance between periodic images (nm)")),
            ("-dz", Option(float, 1, 0, "Z distance between periodic images (nm)")),
            ("-x", Option(vector, 1, 0, "X dimension or first lattice vector of system (nm)")),
            ("-y", Option(vector, 1, 0, "Y dimension or first lattice vector of system (nm)")),
            ("-z", Option(vector, 1, 0, "Z dimension or first lattice vector of system (nm)")),
            #    ("-box",    Option(readBox,     1,        None, "Box in GRO (3 or 9 floats) or PDB (6 floats) format, comma separated")),
            #    ("-n",      Option(str,         1,        None, "Index file --- TO BE IMPLEMENTED")),
            """
        Membrane/lipid related options.  
        The options -l and -u can be given multiple times. Option -u can be
        used to set the lipid type and abundance for the upper leaflet. Option
        -l sets the type and abundance for the lower leaflet if option -u is
        also given, or for both leaflets if option -u is not given. The
        meaning of the number depends on whether option -d is used to set up
        PBC
        """,
            ("-l", Option(self.lipL.append, 1, None, "Lipid type and relative abundance (NAME[:#])")),
            ("-u", Option(self.lipU.append, 1, None, "Lipid type and relative abundance (NAME[:#])")),
            ("-a", Option(float, 1, 0.60, "Area per lipid (nm*nm)")),
            ("-asym", Option(int, 1, None, "Membrane asymmetry (number of lipids)")),
            ("-hole", Option(float, 1, None, "Make a hole in the membrane with specified radius")),
            ("-rand", Option(float, 1, 0.1, "Random kick size (maximum atom displacement)")),
            ("-r", Option(float, 1, 0, "radius of the Membrane sphere")),
            ("-bd", Option(float, 1, 0.3, "Bead distance unit for scaling z-coordinates (nm)")),
            """
        Solvent related options.
        """,
            ("-sol", Option(self.solv.append, 1, None, "Solvent type and relative abundance (NAME[:#])")),
            ("-gas", Option(self.gas.append, 1, None, "Solvent inside type and relative abundance (NAME[:#])")),
            ("-sold", Option(float, 1, 0.5, "Solvent diameter")),
            ("-gasd", Option(float, 1, None, "Gas diameter")),
            ("-solr", Option(float, 1, 0.1, "Solvent random kick")),
            ("-gasr", Option(float, 1, 0.1, "Gas random kick")),
            ("-gden", Option(float, 1, None, "Gas density")),
            ("-excl", Option(float, 1, 1.3, "Exclusion range (nm) for solvent addition relative to membrane surface")),
            """
        Salt related options.
        """,
            ("-salt", Option(str, 1, None, "Salt concentration")),
            ("-charge", Option(str, 1, "auto", "Charge of system. Set to auto to infer from residue names")),
            """
        Define additional lipid types (same format as in lipid-martini-itp-v01.py)
        """,
            ("-alname", Option(self.usrmols.append, 1, None, "Additional lipid name, x4 letter")),
            ("-alhead", Option(self.usrheads.append, 1, None, "Additional lipid head specification string")),
            ("-allink", Option(self.usrlinks.append, 1, None, "Additional lipid linker specification string")),
            ("-altail", Option(self.usrtails.append, 1, None, "Additional lipid tail specification string")),
        ]

    def get_args(self):
        
        # args = sys.argv[1:]
        # output_tk 为字典
        output_tk = InterfaceTK().output_all
        args = []

        for key,value in output_tk.items():
            if key not in ['-u','-sol','-gas','-f']:
                args.extend([key, value])
            else:
                if key == '-u':
                    _extend = []
                    for name,number in value.items():
                        _extend.extend([key,name+":"+number])
                    args.extend(_extend)
                elif key == '-f':
                    pass
                else:
                    value = list(value.keys())[0]
                    args.extend([key,value])
        args.extend(['-z','20'])
        # print(args)


        # if '-h' in args or '--help' in args:
        #     print("\n", __file__)
        #     print(self.desc) or ("\nSomeone ought to write a description for this script...\n")
        #     for thing in self.options:
        #         if type(thing) != str:
        #             print("%10s  %s" % (thing[0], thing[1].description))
        #         else:
        #             print(thing)
        #     sys.exit()
    
        # Convert the option list to a dictionary, discarding all comments
        self.options = dict([i for i in self.options if not type(i) == str])
    
        # Process the command line
        while args:
            ar = args.pop(0)
            self.options[ar].setvalue([args.pop(0) for i in range(self.options[ar].num)])
    
        # Read in the structures (if any)
        # self.tm: list：保存pdb文件的路径
        self.tm = [Structure(i) for i in self.tm]
    
        # absoluteNumbers = not options["-d"]
    
        # 将变量以字典的形式返回, eg:{"self.tm":self.tm}
        variables = [self.tm,self.lipL,self.lipU,self.solv,self.gas,self.usrmols,self.usrheads,self.usrlinks,self.usrtails,self.usrLipHeadMapp,self.usrIndexToLetter]
        args_name = ["tm", "lipL", "lipU", "solv", "gas", "usrmols", "usrheads", "usrlinks", "usrtails", "usrLipHeadMapp", "usrIndexToLetter"]
        variables = dict([(name, var) for name, var in zip(args_name, variables)])
    
        return variables,self.options