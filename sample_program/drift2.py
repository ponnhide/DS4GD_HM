import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
if __name__ == "__main__":
    num_list = []
    for i in range(100):
        num = 0
        box = []
        for j in range(100):
            box.append(str(j+1))     
        
        while len(set(box)) > 1:
            new_box = random.sample(box,50) 
            box = new_box * 2
            num += 1
        num_list.append(num)
    print(np.mean(num_list),np.std(num_list,ddof=1))
    sns.violinplot(data=num_list)
    plt.savefig("violin.pdf")  
