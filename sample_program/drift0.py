import random
if __name__ == "__main__":
    num = 0 
    box = ["red"] * 50 + ["white"] * 50 #これで、"red"が50個、"white"が50個含まれた玉がつくられる
    while len(set(box)) > 1: #box中の要素の種類が2つ以上あったらループし続ける。
        new_box = random.sample(box,50)  #random.sample(list,n)でlistからn個の要素を重複無しにサンプリングができる。
        box = new_box * 2
        num += 1
        print(box) 
    print(num,box)
