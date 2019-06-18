import sys
import matplotlib.pyplot as plt
if __name__ == "__main__":
    positions = [] 
    values    = [] 
    cum_value = 0
    with open(sys.argv[1],"r") as f:
        f.readline()#一行目はseqidとwindow sizeを示すメタ情報なのでスキップ
        for line in f:
            row = line.rstrip().split(",")
            positions.append(int(row[0]))#ポジションは整数に変換 
            cum_value += float(row[1]) 
            values.append(cum_value)
    fig = plt.figure(figsize=(6,2))
    ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    ax.plot(positions,values)
    fig.savefig(sys.argv[1].replace(".txt",".pdf"),bbox_inches="tight")
