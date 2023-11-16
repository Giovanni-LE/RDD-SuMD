mdpparams=dict()
with open("mdp/02-em.mdp",'r') as f:
    for lines in f:
        if lines[0]!=';' and lines[0]!=' ' and lines[0].isalpha():
            mdpparam=lines.strip().replace(" ","").split("=")
            print(mdpparam)
            mdpparams[mdpparam[0]]=mdpparam[1]
print(mdpparams)
mdpparams['include']="ciao"

print(mdpparams)
