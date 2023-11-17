class GmxMdp:
    import sys
    def __init__(self) -> None:
        self.mdp_params=dict()

    def read_mdp(self,mdpfile_in=False):
    #This function read an mdp file
        if mdpfile_in:
            "Give a correct mdpfilename"
            sys.exit(-1)
        with open(mdpfile_in,'r') as f:
            #cicle in the mdp file with trivial check for the correctness of the file
            #check if a line starts with a space
            for lines in f:
                lines=lines.replace.replace(" ","")
                if lines[0]!=';' and lines[0].isalpha():
                    mdpparam=lines.strip().split("=")
                    self.mdp_params[mdpparam[0]]=mdpparam[1]
    def write_mdp(self,mdpfile_out=False):
        if mdpfile_out:
            "Give a correct mdpfilename"
            sys.exit(-1)
        with open(mdpfile_out,'w') as f:
            for key in self.mdp_params.keys():
                f.write(f"{key} = {self.mdp_params[key]}\n")
    def list_mdp_params(self):
        for key in self.mdp_params:
            print(f"{key} = {self.mdp_params[key]}\n")