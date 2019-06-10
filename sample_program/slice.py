hoge = "hogefugahogera"
sliced_hoge = []
size = 2
for i in range(0,len(hoge),size):
    sliced_hoge.append(hoge[i:i+size])
print(sliced_hoge)
