#drift1.py
import random
if __name__ == "__main__":
    num = 0 
    box = []
    for i in range(100):
        box.append(str(i+1)) #boxに1番から100番までの玉を追加。
    while len(set(box)) > 1: #box中の要素の種類が2つ以上あったらループし続ける。
        new_box = random.sample(box,50)  #random.sample(list,n)でlistからn個の要素を重複無しにサンプリングができる。
        box = new_box * 2
        num += 1
    print(num,box)
