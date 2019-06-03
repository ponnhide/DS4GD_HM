import random

if __name__ == "__main__":
    num_list = [] 
    for i in range(1000):
        num = 0
        box = ["red"] * 50 + ["white"] * 50
        while box.count("red") != 100 and box.count("white") != 100:
            new_box = random.sample(box,50)
            box = new_box * 2
            num += 1
        num_list.append(num)
    print(sum(num_list)/1000) 
