import numpy as np

def var(values):
    mean = sum(values)/len(values)
    var  = 0
    for value in values:
        var += (value - mean) ** 2
    var = var / len(values) 
    return var

hoge = [1,2,3,4,5]
fuga = [1,3,5,7,9] 
print(var(hoge)) 
print(var(fuga)) 
