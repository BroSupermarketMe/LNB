# Parse a string for a lipid as given on the command line (LIPID[:NUMBER])
# 解析命令行上给定的脂质字符串 （LIPID[：NUMBER]）
def parse_mol(x):
    '''
    :param x: str: LIPID[:NUMBER]
    :return: (LIPID: str, NUMBER: int or float)
    '''
    l = x.split(":")
    return l[0], len(l) == 1 and 1 or float(l[1])

# Very simple option class 用于处理命令行参数
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self):
        return self.value != None
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]

