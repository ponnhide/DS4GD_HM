hoge = [1,2,3,4,5]
fuga = [1,3,5,7,9] 

hoge_mean = sum(hoge)/len(hoge) 
hoge_var = 0
for value in hoge:
    hoge_var += (value - hoge_mean) ** 2
hoge_var = hoge_var / len(hoge)

fuga_mean = sum(fuga)/len(fuga) 
fuga_var = 0
for value in fuga:
    fuga_var += (value - fuga_mean) ** 2
fuga_var = fuga_var / len(fuga)

print(hoge_var,fuga_var) 
