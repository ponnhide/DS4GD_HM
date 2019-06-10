import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
if __name__ == "__main__":
    population = 100
    num_list   = [] #ループを抜けるまでにかかった回数を保存していくリスト。
    for i in range(100): #100種類の玉が1種類に収束するまでの処理を100回試行する。
        num = 0
        box = []
        for j in range(population):
            box.append(str(j+1))
        while len(set(box)) > 1:
            new_box = random.sample(box,50)
            box = new_box * 2
            num += 1
        num_list.append(num)
    print(np.mean(num_list),np.std(num_list,ddof=1)) #np.meanでlist中の要素の平均値が、np.stdで標準偏差が求められる。ddof=1の意味は、、、
    sns.violinplot(data=num_list)
    plt.savefig("violin.pdf")
