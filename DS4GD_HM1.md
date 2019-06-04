# Python入門

ほんとであれば、Pythonだけで５回ぐらいは授業が欲しいのに、なんと全部で授業は3回しかない。しかもこの授業は「生命動態とデータサイエンス」であって「Python」ではない。なので、Pythonの基本ついては、この１回だけでさらっといきます。

## 特許のMacでPython3をつかう
まず、ターミナルを開いて、下記のコマンドを実行。
````shell
$ bash
$ export PATH=/opt/homebrew/opt/python/libexec/bin:$PATH
````

続けてPythonを実行。
````shell
$ python
Python 3.7.2 (default, Mar 11 2019, 11:23:32)
[Clang 10.0.0 (clang-1000.11.45.5)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>>
````
となればOK。
次は授業で使うモジュール達のinstall

````shell
$ pip install --user biopython
$ pip install --user seaborn
````

特にErrorが出なければOK。

## Introduction of Python
Pythonというのは、ゆとり世代にもわかるように作られたとっても簡単な言語である。基本的にC言語のような古くからあるプログラミング言語っていうのは、低級言語なんて呼ばれたりするのだが、コンピュータにとって敷居が低いように作られていて、ヒトには扱いにくいものが多い。一方でPythonに代表されるスクリプト言語のような高級言語っていうのはコンピュータにとっては敷居が高くて、プログラムを処理、実行するのには時間がかかるけど、ヒトにはとっても理解しやすいように作られた言語である。

現在、Pythnonは生命科学というか科学界を席巻しているプログラミング言語である。日本では日本人が作ったRubyというプログラミング言語があり、日本語のドキュメントが豊富というしょうもない理由でRubyにプログラミング人口がながれたりして、Pythonが流行るのが遅れたけど、そんな日本でも今ではPythonが主流である(Rubyをdisってるわけではなく、RubyはRubyで素晴らしい言語なのだ。）何故こんなに流行ったのか。いろんな理由はあるかもしれないが、大きな理由として可読性の高さと書式の自由度の低さがあるように思う。

Bioinformatics が流行り始めたころ、まずBioinformaticianたちが飛びついたのがPerlと呼ばれるスクリプト言語である。詳しくかかないが、Perlは初期の段階からかなり成熟したスクリプト言語で、テキスト処理にも優れていたため、ゲノムを解析をするバイオインフォマティシャン達の間でめちゃくちゃ流行った。「みんな違ってみんな良い」てきな考えが根底にいったこともあって、初期のPerlプログラマーの間ではどれだけ短くプログラムを書けるかなんて競争も流行ったりした。しかし、「みんな違ってみんな良い」で良かったのは、Bioinformatics という分野がすごーくニッチな分野で、基本的にプログラミングできる奴しかその分野にはいなかったからである。だからこそ、仮にみんな違う書き方でもお互いその内容を理解できたので問題はなかった。しかし、Bioinformaticsという分野が成熟してくると、プログラミングができない人たちもInformatics分野に口を出すようになってきたり、自分たちの手で少しは解析したいとか言い始めるようになってくる。
そうなってくると、「みんな違ってみんな良い」とはならない。むしろ、誰でも同じように書くことができて、誰でもある程度理解できる必要がある。そういう背景で、Perlに取って代わって生命科学を席巻したのがPythonだったということなんだろうなぁとは思ってる。



## Training of Python 

### Pythonの実行方法
ターミナルに```python```と打ち込んでやれば、
````shell
$ python
Python 3.7.2 (default, Mar 11 2019, 11:23:32)
[Clang 10.0.0 (clang-1000.11.45.5)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>>
````
こんな感じで、インタプリタと呼ばれるものが立ち上がる。インタプリタは簡単な計算や動作の確かめには便利。ただ、大きなプログラムを動かすにはちと不便。そう言った場合にはファイルにプログラム（スクリプト）を書き込み、そのファイルを実行することになる。次ページ参照。

<div style="page-break-after: always;"></div>

### Hello World 

テキストエディタを立ち上げて（間違えてもWordとか立ち上げないでね。Wordにプログラム書き始めるのだけはマジでやめてね。）以下のプログラムを書き込こもう。
````python:hello.py
print("Hello world")
````
そしたら"hello.py"で保存。そのまま、Terminal上で以下のコマンドを実行。(ちゃんとhello.pyを作成したディレクトリまで移動してね。)
````shell
$python hello.py
````
"Hello world"って表示されたかな？。こんなの簡単じゃないかと思うかもしれんが、未だにネットとか調べると
````python
print "Hello world"
````
とか書いてる古典が散乱している。この書き方はPython2系の書き方である。今更こんな書き方を覚える必要はないというかPython3ではエラーを吐かれて動かない。ところで、なぜ2系から3系にするにあたり変更したのか。書くと長くなるので気になった人は調べてみよう。

また、"Python Hello world"でググると、下記のようなプログラムを見かけることがある。
`````Python
if __name__ == "__main__":
	print("Hello world")
`````
これは正しいのだが、```if __name__ == "__main__"```の意味にについて話すと長くなるので授業では省きます。基本的に、Pyhtonプログラムを書く際には```if __name__ == "__main__":```をつけた方がいい。ただ、今はどっちでもおkです。



### 変数、代入
Pythonには驚くべきことに変数宣言も型宣言が存在しない。なので、代入プロセス自体が変数宣言とほぼ同義である。で、先週までRを使っていた人たちに気をつけて欲しいのは、代入演算子が "<-" でなく "=" であるということである。実を言うとRも"="はつかえるけどね。"="は代入である。これさえ理解すればプログラミングはもう怖くない。
````Python
a = 1
b = 2
b = a
a = 3
print(a) 
print(b)
````
さて、上記のa, bの値は何になるでしょう。正解は print(a) が 3で print(b)が2。簡単な問題にみえるけど実は奥が深い。



### Pythonを使った四則演算
以下のコードをPython3で実行してみよう。
````Python
a = 5
b = 2
print(a + b) 
print(a - b)　
print(a / b)
print(a * b)
print(a % b)
print(a // b) 
print(a ** b)
````
各々がどんな意味を持つかわかったかな？。実をいうとPython2で同じコードを実行すると結果が変わる。気になった人はPython2で実行してみよう。

<div style="page-break-after: always;"></div>

### List型
Listは複数のデータを格納するための型。C言語とかと違って格納できる要素の型に制限はなく、異なる型の要素を同じlistに格納することもできる。例えば、下記のように数値と文字列を同じlistにいれることも可能
````Python
hoge = [1,"hoge"]
````
List中の要素の取り出し方、追加、変更は以下のように行う。
````Python
hoge = [1,2,3]
print(hoge[0]) #1つ目の要素
print(hoge[1]) #2つ目の要素
hoge.append(4) #4つ目の要素の追加
hoge[2] = 5    #3つめの要素の変更
print(hoge)
````
Rと違ってindexが0から始まる点に注意。n番目の要素を取り出したい場合にはn-1を用いる。一般的にプログラミングの世界では0から数える方が一般的。Rとかの方が珍しい。



**※List型のコピー**

List型を扱う上で気をつけないといけないのは以下のような例。

`````Python
hoge = [1,2,3]
fuga = hoge
fuga[0] = 4
print(hoge)
print(fuga)
`````
さて、なんて表示されるかな。ここを乗り切れれば、もうPythonなんて楽勝である。
基本的にlist型のコピーを行う際には、
`````Python
hoge = [1,2,3]
fuga = hoge[:]
fuga[0] = 4
print(hoge)
print(fuga)
`````
とすればOK。ただし多重リスト(list中の要素についても別のlistなっているlist)の場合はこの限りではない。多重リストをコピーしたい場合には「**辞書型のコピー**」の項目を参考にしてほしい。



### 条件分岐と繰り返し処理
みんな情報基礎で理解していると信じているので、細かい部分は省略。条件分岐と繰り返し処理は組み合わせて使うことが多い。例えば、各genomeのGC含量の格納したlistから、GC含量が0.6以上のゲノムの数を数えようと思ったらこんな感じ。
````Python
GC_ratios = [0.52,0.48,0.47,0.63,0.71,0.50,0.42] #example data
num = 0
for gc_ratio in GC_ratios:
    if gc_ratio >= 0.6:
    	num += 1
print(num) 
````

何番目の要素がGC含量≥0.6か知りたかったら、以下のようにすれば、要素のindexと値を同時に取得できる。
````Python
for i in range(len(GC_ratios)):
    if gc_ratio >= 0.6:
    	print(i,GC_ratios[i])
print(num) 
````
以下のように書いてもおk。
````Python
for i gc_ratio in enumerate(GC_ratios):
    if gc_ratio >= 0.6:
    	print(i,gc_ratio)
print(num) 
````

ここまでfor文を使った繰り返しだったけど、ある条件を満たしている間は同様の処理を繰り返すみたいなときは while 文を使う。例えば、こんな感じ。
````Python
import random
num  = 0
rand = 0 
while rand < 0.5:
	rand = random.random()
	num += 1
print(num,rand) 
````
ちょっと、新しいことを出してしまったので上記のプログラムを解説します。
上記のプログラムでは、while文を使って変数randが0.5未満の場合は乱数を発生させ続けるという処理をしている。ただし、実はPythonそのものには乱数を発生させる機能はない。そこで、最初の行の```import random```で乱数を発生させるためのプログラム（モジュール）を外部から呼び出している。これが「モジュールの読み込み」とか言われるやつ。こうすることで、random.〜とやれば乱数生成に関連する様々なメソッドが使えるようになる(〜にはメソッドの名前が入る）。上記のプログラムではradnomモジュールにあるrandom メソッド```random.radnom()```を使って0~1の乱数を発生させている。例えば他にも、```random.radint(0,10)```を使えば、0〜10までの整数をランダムに発生させることができる。**こんな風にPythonでは様々なモジュールを読み込んで便利な機能を使うことができるのだ！**



### 辞書型
keyとvalueが組みなるデータ型。異なる情報を紐づけたいときに利用する。例えば、
````Python
first_last = {"Hideto":"Mori","Haruo":"Suzuki"}
print(first_last["Hideto"]) 
````
みたいな使い方。keyを元にvalue("Hideto"）アクセスできたことがわかる。valueの値を変えたり、新たなkeyを追加することも可能。
````Python
first_last["Hideto"] = "Yamamoto"
print(first_last["Hideto"]) 
first_last["Masaru"] = "Tomita"
print(first_last)
````
Hidetoに紐付いたValueがYamamotoになっていること、新たにMasaruのkeyが追加されたことを確認しよう。



**辞書型のコピー**

辞書型のコピーについてもListのコピーと同じことが起こるので注意。辞書や多重Listを他の変数にコピーしたい場合にはcopyモジュールを使うと良い。

````Python
import copy
first_last = {"Hideto":"Mori","Haruo":"Suzuki"}
hoge = first_last    #=をつかったコピー（代入）
fuga = copy.deepcopy(first_last)    #copy.deepcopy()を使ったコピー
hoge["Hideto"] = "Yamamoto"
fuga["Haruo"]  = "Sato"

print(first_last)
print(hoge) 
print(fuga) 
````
どんな結果になったかな？。確かめてみよう。



### Set型
あまり使わないけど覚えておくと便利なのがSet型。list中に含まれる要素の種類および種類数を調べるときとかに使う。例えば、
````Python
hoge = [1,2,1,3,3]
fuga = set(hoge)
print(fuga) 
````
とやると、```set([1, 2, 3])```が得られる。こんな風にListの重複をなくして、中にある要素の種類だけ取り出してくれる。

<div style="page-break-after: always;"></div>

### 生命動態のシミュレーションと簡単な統計処理

さて、ここからやっと実践的な話。せっかくPythonの基本について理解したので、簡単な生命動態のシミュレーションについて触れようと思う。
以下の問題について、その結果がどうなるか考えてみてください。

**問題**
50個の赤い玉と50個の白い玉があるとき、

1. 無作為に50個の玉を選ぶ。
2. 選んだ玉を赤と白の比を保ったまま2倍に増やす。(50個->100個)

1と2の動作を延々と繰り返していくと、最終的に100個の玉はどういった状態に落ち着くでしょうか。

上の問題は確率と極限を使えば計算で結論を導くことができるのだが、解説が面倒なのでここでは飛ばします。結論からいうと必ず赤が100個、白が100個の状態に落ち着く。プログラムを作成して実際に確かめてみよう。

````Python
import random
if __name__ == "__main__":
    num = 0 
    box = ["red"] * 50 + ["white"] * 50 これで、"red"が50個、"white"が50個含まれた玉がつくられる
    while len(set(box)) > 1: #box中の要素の種類が2つ以上あったらループし続ける。
        new_box = random.sample(box,50)  #random.sample(list,n)でlistからn個の要素を重複無しにサンプリングができる。
  	    box = new_box * 2
  	    num += 1
    print(num,box)
````
さて、上記のプログラムはboxのlistの中身が全て赤玉 or 白玉になるまで1と2の試行を繰り返すものである。何度実行してもらっても構わないが、何度実行しようと必ずwhileループを抜けて白玉100個もしくは赤玉100個の状態になる。



問題をちょっと変えて、1番から100番までの番号がついた100個の玉を対象に、先と同じ処理を行ってみよう。
````Python
import random
if __name__ == "__main__":
    num = 0 
    box = []
    for i in range(100):
        box.append(str(i+1)) boxに1番から100番までの玉を追加。
    while len(set(box)) > 1: #box中の要素の種類が2つ以上あったらループし続ける。
        new_box = random.sample(box,50)  #random.sample(list,n)でlistからn個の要素を重複無しにサンプリングができる。
  	    box = new_box * 2
  	    num += 1
    print(num,box)
````

これでも、やっぱり必ず100個の玉は1つの数字に収束する。何度やろうと収束する。
せっかくなのでこのプログラムを100回繰り返して、1つの数字に収束するまでかかるループの回数の平均と分散を求めてみよう。
以下のプログラムを実行すると、100試行分の平均、分散およびループ回数の分布を示したviolin plot を得ることができる。

````Python 
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
    print(np.mean(num_list),np.std(num_list,ddof=1))　#np.meanでlist中の要素の平均値が、np.stdで標準偏差が求められる。ddof=1の意味は、、、
    sns.violinplot(data=num_list)
    plt.savefig("violin.pdf")
````

で、結果はどうでも良いとして、これのどこが**生命動態のシミュレーション**なのか。それは授業で話します。



### 課題
最後のプログラムにおける、population = 100の部分を200,500,1000に変えてプログラムを実行してみよう。各々の場合で出力されるpdfをwordに貼り付けて、結果に対する考察と授業の感想も書いて提出してください。



## おまけ
###関数
プログラミングをやる上で絶対に知っておいた方がいいものの1つが関数である。関数とは特定の処理に名前をつけてまとめたものである。例えば、numpyを使わずに2つのList各々の分散を求めるためのプログラムを考えてみよう。
````Python 
import numpy as np
hoge = [1,2,3,4,5]
fuga = [6,7,8,9,10] 

hoge_mean = sum(hoge)/len(hoge) 
hoge_var = 0
for value in hoge:
	hoge_var += (value - hoge_mean) ** 2
hoge_var = hoge_var / len(hoge)

fuga_mean = sum(fuga)/len(fuga) 
fuga_var = 0
for value in hoge:
	fuga_var += (value - fuga_mean) ** 2
fuga_var = fuga_var / len(fuga)

print(hoge_var,fuga_var) 
````
プログラムをみたら気づいたと思うけど、fuga_mean からはhogeに対して処理をfugaに対しても繰り返しているだけである。これは、無駄そうだし変数名も増えていってプログラムが汚くなってしまう。ということで、こういう時は関数を使えばよい。
````Python 
import numpy as np
def var(values):
	mean = sum(values)/len(values)
	var  = 0
	for value in values:
		var = += (value - mean) ** 2
	var = var / len(values) 
	return var
hoge = [1,2,3,4,5]
fuga = [6,7,8,9,10] 
print(var(hoge)) 
print(var(fuga)) 
````
さて、関数の役割はわかったかな？










