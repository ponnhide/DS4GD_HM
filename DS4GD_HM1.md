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
とか書いてる古典が散乱している。この書き方はPython2系の書き方である。今更こんな書き方を覚える必要はないというかPython3ではエラーを吐かれて動かない。ところで、なぜ2系から3系にするにあたりこのような変更をしたのか。書くと長くなるので気になった人は調べてみよう。

また、"Python Hello world"でググると、下記のようなプログラムを見かけることがある。
`````Python
if __name__ == "__main__":
	print("Hello world")
`````
これは正しいのだが、```if __name__ == "__main__"```の意味にについて話すと長くなるので後回し。基本的に、Pyhtonスクリプトを書く際には```if __name__ == "__main__":```をつけた方がいい。ただ、今はどっちでもおk。



### 変数、代入
Pythonには驚くべきことに変数宣言も型宣言が存在しない。なので、最初の代入プロセス自体が変数宣言のようなものである。で、先週までRを使っていた人たちに気をつけて欲しいのは、代入演算子が "<-" でなく "=" であるということである。実を言うとRも"="はつかえるけどね
**"="は代入である。これさえ理解すればプログラミングはもう怖くない。**

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
List中の要素の取り出し方、追加、変更、要素数の確認は以下のように行う。
````Python
hoge = [1,2,3]
print(hoge[0]) #1つ目の要素
print(hoge[1]) #2つ目の要素
hoge.append(4) #4つ目の要素の追加
hoge[2] = 5    #3つめの要素の変更
print(len(hoge)) #len（）を使うことでリストの要素数がわかる。
print(hoge)
````
Rと違ってindexが0から始まる点に注意。n番目の要素を取り出したい場合にはn-1を用いる。一般的にプログラミングの世界では0から数える方が一般的。Rとかの方が珍しい。これ以降は、0番目と言った場合には1つめの要素。n番目の要素といった場合にはn+1つ目の要素のこととする。ついでにリスト同士を連結する場合は+が使える。ただ、```.extend()```を使う方が一般的。またあまり知られていないが、\*を使うとで同じ要素を繰り返すリストが作れる。
````
hoge=[1,2,3]
fuga=[4,5,6]
hoge = hoge + fuga
print(hoge) 
hoge.extend(fuga)
print(hoge)
fuga = fuga * 3 
print(fuga) 
````


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

実行してくれたら分かると思うが、fugaの中身を変えたはずなのに、hogeの中身まで変わってしまう。この理由については長くなるので省くが、これでは困るので、list型のコピーを行う際には、

`````Python
hoge = [1,2,3]
fuga = hoge[:]
fuga[0] = 4
print(hoge)
print(fuga)
`````
とすればOK。ただし多重リスト(list中の要素も別のlistなっているlist、行列とか)の場合はこの限りではない。多重リストをコピーしたい場合には「**辞書型のコピー**」の項目を参考にしてほしい。どうしてもここら辺の理由をちゃんと知りたいひとは、「Pythonの参照の値渡し」でググってみよう。あと、自分のQiitaの宣伝になってまうけど、なんとなく分かればいいかなって人は<https://qiita.com/ponnhide/items/cda0f3f7ac88262eb31e>を読むといいかも。



###文字列型
そういえば、初回の授業で文字列型の話忘れた。基本的にどんな文字も数字もダブルコーテーションもしくはシングルコーテーションで囲ってやれば文字列型扱いになる。文字列はリストとほぼ同じように扱えるのだが、一点だけ違うのは要素の変更ができないということである。よって以下のようなことはできない。
````Python
hoge="HIDETO"
hoge[1] = "O"
````
上記のようなことをしたい場合は、一旦list型に変更して、最後に文字列型に戻すのが楽。
````Python
hoge = list(hoge) 
hoge[1] = hoge
hoge = str(hoge) 
print(hoge)
````
文字列の連結には、listと同じように+記号をつかう方法と```.join()```を使う方法がある。

````Python
hoge = "MORI"
fuga = "HIDETO"
hogera = hoge + fuga
fugara = " ".join([hoge,fuga]) 
print(hogera)
print(fugara)
````
join()の使いかたは分かっただろうか。hogeとfugaがスペースで区切られて連結されたことが確認できたと思う。スペース以外にも、```",".join()```とすればカンマで区切りに、```":".join()```とすればコロン区切りになる。```.join()```はよく使うので覚えておくと便利。



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
何となく、動作はわかったと思うけどちっと解説。上の例ではfor文を使うことで変数gc_ratio にGC_ratios中の各値が1つづつ代入されて、for文の中で処理されている。こんな風にfor文をつかうと、リスト中の値を１つづつ取り出して処理を行うことができる。え？リストに限らずイテレートできる型なら、全部for文で処理できるって？。まぁうん、でもイテレーターとか言っても難しいし、ややこしいからな。うん。

単に順番にリストの中身を取り出すんじゃなくて、取り出した要素が「何番目の要素であるという情報」も知りたい場合、以下のようにすればいい。
````Python
for i in range(len(GC_ratios)):
    if GC_ratios[i] >= 0.6:
    	print(i+1,GC_ratios[i])
print(num) 
````
```range(n)``` (nは正の整数）としてやると、0からn-1まで数字をループ1回ごとに順番に生成(ジェネレイト）してくれるジェネレーターというものが作られる。したがって、上の場合```i```には0から6までの数字が順番に代入されるので、```GC_ratios[i]```とすれば、```i```番目の要素を取り出すことができる。

以下のように書いても同じことができる。こっちの書き方をするとチョットデキルひとに見えるかも！
````Python
for i gc_ratio in enumerate(GC_ratios):
    if gc_ratio >= 0.6:
    	print(i,gc_ratio)
print(num) 
````
解説が面倒になってきたので、"enumerate python"でググってくれ。うん。

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
上記のプログラムでは、while文を使って変数randが0.5未満の場合は乱数を発生させ続けるという処理をしている。ただし、実はPythonそのものには乱数を発生させる機能はない。そこで、最初の行の```import random```で乱数を発生させるためのプログラム（モジュール）を外部から呼び出している。これが「モジュールの読み込み」とか言われるやつ。こうすることで、random.〜とやれば乱数生成に関連する様々な関数が使えるようになる(〜には関数名が入る）。上記のプログラムではradnomモジュールにあるrandom 関数```random.radnom()```を使って0~1の乱数を発生させている。例えば他にも、```random.radint(0,10)```を使えば、0〜10までの整数をランダムに発生させることができる。**こんな風にPythonでは様々なモジュールを読み込んで便利な機能を使うことができるのだ！**



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
```.keys()```や、```.values()```を使ってkeyだけ、valueだけを取り出すこともできる。ただし```.keys()```や、```.values()```は```range()```と同様に順番にしたがってkeyとvalueを1つづつ返すだけなので、中身全体を確認したい場合は```list()```で囲んでやる必要がある。
````Python
print(list(first_last.keys()))
print(list(first_last.values())
````
あと、古典みたいなwebサイトにはPythonの辞書はkeyを追加したときに、順序が保存されないとか書いてあることがあるが、そんなことはない。Python3.6以降を使っていればちゃんと保存される。


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
とやると、```set([1, 2, 3])```が得られる。こんな風にListの重複をなくして、中にある要素の種類だけを教えてくれる。
<div style="page-break-after: always;"></div>


