## 生命動態(遺伝的浮動）のシミュレーションと簡単な統計処理
さて、ここからやっと実践的な話。せっかくPythonの基本について理解したので、簡単な生命動態のシミュレーションについて触れようと思う。以下の問題について、その結果がどうなるか考えてみよう。

**問題**
50個の赤い玉と50個の白い玉があるとき、

1. 無作為に50個の玉を選ぶ。
2. 選んだ玉を赤と白の比を保ったまま2倍に増やす。(50個->100個)

1と2の動作を延々と繰り返していくと、最終的に100個の玉はどういった状態に落ち着くでしょうか。

上の問題は確率と極限を使えば計算で結論を導くことができるのだが、解説が面倒なのでここでは飛ばします。結論からいうと必ず赤が100個、白が100個の状態に落ち着く。プログラムを作成して実際に確かめてみよう。

````Python
#drift0.py
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
````
さて、上記のプログラムはboxのlistの中身が全て赤玉 or 白玉になるまで1と2の試行を繰り返すものである。何度実行してもらっても構わないが、何度実行しようと必ずwhileループを抜けて白玉100個もしくは赤玉100個の状態になる。

問題をちょっと変えて、1番から100番までの番号がついた100個の玉を対象に、先と同じ処理を行ってみよう。
````Python
#drift1.py
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
せっかくなのでこのプログラムを100回繰り返して、1つの数字に収束するまでかかるループの回数の平均と分散を求めてみよう。以下のプログラムを実行すると、100試行分のループ回数の平均、分散およびその分布を示したviolin plot を得ることができる。

````Python 
#drift2.py
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
````
で、結果はどうでも良いとして、これのどこが**生命動態のシミュレーション**なのか。実はこれ、[**遺伝的浮動**](https://ja.wikipedia.org/wiki/%E9%81%BA%E4%BC%9D%E7%9A%84%E6%B5%AE%E5%8B%95)のシミュレーションなんですね。伝えたかったことは、選択圧がかかっていなくとも世代数をかければ必ず集団は特定の1個体からの子孫によって固定されるという事である。つまり、集団における遺伝的多様性が失われることを意味する。特にこの遺伝的浮動の影響は、集団のサイズが小さいほど大きくなる。また、遺伝的浮動の影響は選択が掛かっている場合においても存在する。授業では触れないが、できる人は上記プログラムを"1"番の玉が他の玉より2倍選ばれやすくなるように書き換えて、再実行してみよう。"1"番に収束する確率は大きくなるものの、他の数字にも収束することが依然としてあることが確認できるはずだ。

こんな風にごく簡単なプログラムで、遺伝的浮動の様子をシミュレーションをすることができる。にも関わらず、進化において遺伝的浮動の影響が議論されるようになったのはたかだが40-50年前のことであり、それまで生命の設計は自然選択によって決定されるという考えが主流だった。しかし、実際には現生の生命の設計は偶然によって固定された要素を多分に含んでいるのである。

**授業内課題**
drift_2.pyにおける、population = 100の部分を50, 200, 500に変えて実行してみよう。何か法則に気づいたかな。



## 関数（授業では飛ばします）
プログラミングをやる上で絶対に知っておいた方がいいものの1つが関数である。関数とは特定の処理に名前をつけてまとめたものである。例えば、numpyを使わずに2つのList各々の分散を求めるためのプログラムを考えてみよう。え、分散の求め方がわからない？。。。。。。。。[wikipedia](https://ja.wikipedia.org/wiki/%E5%88%86%E6%95%A3_(%E7%A2%BA%E7%8E%87%E8%AB%96))

````Python 
#var0.py
hoge = [1,2,3,4,5]
fuga = [1,3,5,7,9] 

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
プログラムをみると気づくと思うが、fuga_mean からはhogeに対して行った処理をfugaに対しても繰り返しているだけである。これは、無駄そうだし変数名も増えていってプログラムが汚くなってしまう。ということで、こういう時は関数を使えばよい。

````Python 
#var1.py
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
````
さて、関数の役割はわかったかな？**基本的に関数の中で定義した変数は関数の外では使えない**ので、 ```return var```で関数内で行った処理の結果である```var```を返している。返って来た結果は、```hoge_var=var(hoge)```のように新たな変数に代入することもできる。

で、関数の話をしたので、先週でてきた```if __name__ == __main__```についてもついでに解説しちゃおう。これまでにも話したようにPythonではモジュールとして別ファイルに書かれた関数を利用することができる。例えば```import random```を行うと、random.pyに書かれている関数を```random.~```として使うことができる。これは、自分が作ったファイルについても同じである。(そもそもrandom.pyはどこにあるのか。実はrandom.pyがあるディレクトリまでPATHが通ってるので、ユーザ側がファイルの場所を明示しなくても読み込むことができる。というか、PATHが通ってない場所のファイルはモジュールとして読み込めない。気になった人はPYTHONPATHでググろうね。)

では、sample_program中に、import_example0.pyとimport_example1.pyというファイルがあるので、両方を実行してみよう。実行結果が違うのが分かったかな？この違いはどこから来るのか。ファイルの中身を見るとわかるが、
````Python
#import_example0.py
import var1
print(var1.var([1,2,3,4,5])) 
````
````Python
#import_example1.py
import var2
print(var2.var([1,2,3,4,5])) 
````
import_example0.pyでは、同じディレクトリ中のvar1.pyを読み込んでいるのに対して、import_example1.pyはvar2.pyを読み込んでいる。では、var1.pyとvar2.pyの違いは何か。それが```if __name__ == __main__```である。
````Python 
#var2.py
import numpy as np

def var(values):
    mean = sum(values)/len(values)
    var  = 0
    for value in values:
        var += (value - mean) ** 2
    var = var / len(values) 
    return var

if __name__ == "__main__":
    hoge = [1,2,3,4,5]
    fuga = [1,3,5,7,9] 
    print(var(hoge)) 
    print(var(fuga)) 
````
さて、どうしてimport_example0.pyとimport_example1.pyの実行結果に違いがあったのか分かったかな。実行結果を見るとimport_example0.pyではvar1.py中の```print(var(hoge))```と```print(var(fuga))```が```print(var_1.var([1,2,3,4,5]))```より前に実行されている 。一方で、 import_example1.pyでは、var2.py中の ```print(var(hoge))```と```print(var(fuga))```は実行されず、```print(var2.var([1,2,3,4,5])) ```しか実行されない。これは、Pythonファイルが実行ファイルとして呼び出された場合には```__name__```が```"__main__"```となるのに対して、モジュールとして呼び出された場合には、```__name__```がモジュール名になるためである。したがって、var2.pyをモジュールとしてimportした場合には、```if __name__ == __main__```以下は条件を満たさないから実行されないけど、```python var2.py```のように実行した場合には```if__name__ == "__main__"```以下まで実行されるのだ！。わかったかな？。



## DNA/Protein sequence analysis wiith Python
ここでは、BioPythonを使ってゲノム配列データを解析する手法について扱う。BioPythonはめちゃくちゃ機能が多い割に、解説が少ない泣。今日の授業では、Fastaファイルの基本的な扱い方とGC含量、GC skewの移動プロットまでを扱う。

### Fastaの読み込み->逆相補鎖配列->GCratio->GCskew

#### with Biopython
````Python 
#sample1_w_bio.py
from Bio import SeqIO 
from Bio import SeqUtils 
if __name__ == "__main__": 
    for fna in SeqIO.parse("../sample_genome/sample1.fasta","fasta"):
        print(fna.id)
        print(fna.seq) 
        print(fna.seq.reverse_complement()) 
        print(SeqUtils.GC(fna.seq)) 
        print(SeqUtils.GC_skew(fna.seq)[0])
````
一応、1つ1つ解説していく。
**from Bio import SeqIO**   
Bio(python)から、様々な形式で書かれるゲノム情報を読み込むためのモジュールSeqIO (Iはinput, Oはoutputの意味）を読み込んでいる。```import Bio.SeqIO``` といった書き方も可能。

**from Bio import SeqIO**   
Bio(python)から、生体分子配列 (DNA/RNA配列）に対して様々な生命科学計算を行ってくれるモジュールSeqUtilsを読み込んいる。```import Bio.SeqUtils``` といった書き方も可能。

**SeqIO.parse("../sample_genome/sample1.fasta","fasta")**  
```SeqIO.parse```を使うとファイル内の配列情報を1配列ごとに取り出すことができる。最初の引数にfileの名前、次にフォーマット(Fasta, genbank, gff等等）を指定する。念の為に説明すると"../"は１つ上のディレクトリを示す相対PATHだぞ。同じようなことをしてくれる関数にSeqIO.readという関数がある。しかし、こちらはファイル内の配列数が1つのときしか使えないので、あまり使い道がない。基本的に、FastaもGenbankも複数の配列の情報が入っていることが多い。ヒトであれば、各染色体の配列情報だったり、バクテリアであればまだゲノムとプラスミド、コンプリートゲノムが分かっていない生物なら各contigの配列等々。

**fna.id**  
Fasta fileのフォーマットは以下のようになっている。ここでいう```HOGE```や```FUGA```にあたる部分が、idである。
````
>HOGE
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTG
GTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGAC
>FUGA
GTGTTGCCAACTCGAAGGCTCTGCTCACCAATGTACATGGCCTTAATCTGGAAAACTGGCAGGAAGAACTGGCGCAAGCC
AAAGAGCCGTTTAATCTCGGGCGCTTAATTCGCCTCGTGAAAGAATATCATCTGCTGAACCCGGTCATTGTTGACTGCAC
````

**fna.seq.reverse_complement()**  
逆相補鎖の配列を作成してくれる。

**SeqUtils.GC(sequence)**  
一本鎖DNA配列中のGC含量を計算してくれる。

**SeqUtils.GC_skew(sequence, windowsize=100)**  
一本鎖DNA配列中のGC skewを計算してくれる。**GC skewは一本鎖DNA配列におけるG含量とC含量のバイアスを表す指標**。(G-C)/(G+C)の式で表される。この値が正であれば、その配列中でGにバイアスがかかっていることになるし、負であればCのバイアスが大きいことになる。お節介なことに、この関数は2つ目の引数(windowsizeの部分）に数字を入れると、そのサイズのwindowでGC_skewの移動平均を取ってくれる（何も入れなかった場合には勝手にwindowsize=100の移動平均になる）。なぜこの機能が```SeqUtils.GC```にはないのか。謎である。

#### without Biopython
ぶっちゃけ上記のプログラムぐらいなら、Biopythonを使わなくたってできる。シンプルなプログラムならBiopyhtonを使わない方が処理も早いので、自分で書けるようになることをお勧めする。

````python
#sample1_wo_bio.py
def read_fasta(fname):
    seq_dict = {} 
    with open(fname,"r") as f:
        for line in f:
            if line[0] == ">":
                seq_id = line.rstrip()[1:]
                seq_dict[seq_id] = ""
            else:
                seq_dict[seq_id] += line.rstrip().upper()#念のため大文字化 
    return seq_dict

def calc_GC(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) 
    
def calc_GC_skew(seq):
    G_num = seq.count("G")
    C_num = seq.count("C")
    return (G_num - C_num)/ (G_num + C_num)
    
if __name__ == "__main__": 
    seq_dict = read_fasta("../sample_genome/sample1.fasta")
    for seq_id in seq_dict.keys(): 
        seq = seq_dict[seq_id]
        print(seq)
        print(seq.translate(str.maketrans("ATGC","TACG"))[::-1]) 
        print(calc_GC(seq))
        print(calc_GC_skew(seq))
````
一部だけ、解説
**read_fasta(filename)**  
Fastaを読み込んで、SeqIDをkey, 配列をvalueに持つ辞書を返してくれる関数を作っている。
>**with open(fname,"r") as f**  
>ここではfileを読み込んでいる。fnameにはファイル名、fには任意の変数名が入る。この場合「ファイル(fname)をfとして読み込みますよー」といったところである。"r"の部分には基本的は"r"か"w"が入る。"r"は読み込みモード、"w"が書き込みモードである。間違ってもすでに存在するファイルを"w"で開いてはいけない。中身が消されて初期化されてまうからな。既にあるファイルに新たな書き込みを加えたい場合には"a"を使おう。
>**for line in f:**  
>こうすることで、ファイルを1ループに1行づつ読み込むことができる。基本的に途中から読み込んだりとかはできない。
>**if line[0] == ">"以下の解説**  
>各行の一文字目が">"だったら、```rstrip()```を使って改行文字を取り除いた後、```[1:]```で">"以降の文字列を```seq_id```に代入し、```seq_dict```に新たなkeyとして加えている。先頭が”>"以外の行は配列情報なので、rstrip()で改行文字を取り除いた文字列を```seq_dict[seq_id]```に随時連結している。

**seq.translate(str.maketrans("ATGC","TACG"))[::-1]**  
これは公式として覚えて欲しい。こうすれば、逆相補鎖になるのだ。なるものはなる。単に相補鎖が欲しい場合は、```seq.translate(str.maketrans("ATGC","TACG"))```でおk。



### GC含量の移動プロット
GC含量の移動プロットを書いてみよう。残念ながら、Biopythonには1行で移動プロットをしてくれるような関数はない。なので、自分で作るしかない。まぁそんな難しいコードでもないので、書いてみよう。

まず移動プロットを書くためには、ゲノムをウインドウサイズごとに分ける必要がある。ようは文字列をn文字ごと区切れれば良いのだが、どのように行うかというと、
````Python
#slice.py
hoge = "hogefugahogera"
sliced_hoge = []
size = 2
for i in range(0,len(hoge),size):
    sliced_hoge.append(hoge[i:i+size])
print(sliced_hoge)
````
まぁこんな感じで、sizeのところを任意の数に変えればsize毎に文字列を分割することができる。では早速使ってみよう。

````Python 
#GC_slide1.py
from Bio import SeqIO 
from Bio import SeqUtils 
if __name__ == "__main__": 
    size = 10000
    for fna in SeqIO.parse("../sample_genome/sample2.fasta","fasta"): 
        results = []
        seq = fna.seq
        for i in range(0,len(seq),size):
            GC_ratio = SeqUtils.GC(seq[i:i+size])
            results.append(GC_ratio)
        print(fna.id,results) 
````
まぁ、これだと単に10,000 bpごとのGC含量のlistが表示されるだけなので、ちゃんと結果をプロットしてみよう。



#### 図のプロット
ここからは個人的な流儀なのだが、どんな結果も基本的に一度テキストファイルに書き出した方が良い。プログラミングに慣れてくるとデータを書き出さず、1つのスクリプトの中でデータの処理から結果の図示まで完結させがちになる（僕自身、面倒なときはこうしてしまう時がある。）ただ、これはあまり良くない。理由は簡単で、他の人が結果を確認できないからだ。世の中の人誰もがプログラムを走らせられる訳ではない。エクセルしか使えない人だってたくさんいるっていうかそっちの方が多い。で、研究をつづけていけばそういう人たちとデータを共有する機会が増えてくる。そういう時は、図だけ送っても相手には生の数値がわからない。一方、エクセルで再現できるように整形した生データがあれば、相手は数値を確認できるし、自分の手で同じ図を再現することだってできる。なので図を作る際は、元となる数値データをテキストに書き出すためのスクリプトと、そのテキストを読み込んで図をプロットするスクリプトは別々に作成するようにしよう。

ということで、先の結果をテキストに書き出してみよう。
````Python 
#GC_slide2.py
from Bio import SeqIO 
from Bio import SeqUtils 
if __name__ == "__main__": 
    size = 10000
    for fna in SeqIO.parse("../sample_genome/sample2.fasta","fasta"): 
        results = []
        seq = fna.seq
        for i in range(0,len(seq),size):
            GC_ratio = SeqUtils.GC(seq[i:i+size])
            results.append(GC_ratio)
        
        with open("{}_GCratio.txt".format(fna.id),"w") as o:
            o.write(fna.id + "," + str(size) + "\n")
            for i in range(len(results)):
                o.write(str(i*size+1) + "," + str(results[i]) + "\n")
````
もう私には無理ぽ。今日は帰って寝よ。って感じの子が出てきそうなので、ちょっと解説。

**with open("{}_GCratio.txt".format(fna.id),"w") as o**  
```sample1_wo_bio.py```のときとは違って、ここではファイルを書き込みモードで開いている。```"{}_GCratio.txt".format(fna.id)```で使われている```format```は{}の中に変数を埋め込むために使われれる。詳しい使い方は"Python format"でググってみよう。

**o.write(fna.id + "," + str(size) + "\n"), o.write(str(i*size+1) + "," + str(results[i]) + "\n")**  
o.writeでファイルに文字を書き込むことができる。この時に気をつけたいのが、改行したい場合には最後に"\n"をつけること。これは改行コードと呼ばれる文字で、これをつけないと改行が行われないので注意。

さて上記のプログラムを動かすと、
````
sample,10000
1,52.07
10001,49.94
20001,52.62
30001,53.23
40001,52.77
50001,51.56
60001,55.58
70001,53.6
````
のようなファイルができたはず。次は、このファイルを読み込んで結果をプロットしてみる。
Pythonで図を作る際には、matplotlibというモジュールを使う。matplotlibにできないプロットはないと言ってもいいほど自由度の高いライブラリである。しかし、その自由度の高さはもはや悪魔である。使いこなしたい人はとりあえず、
https://qiita.com/skotaro/items/08dc0b8c5704c94eafb9
を読んでみよう。ただ、今回は単にGC含量をプロットしたいだけなので適当にいく。

````Python
#plot.py
import sys
import matplotlib.pyplot as plt
if __name__ == "__main__":
    positions = [] 
    values = [] 
    with open(sys.argv[1],"r") as f:
        f.readline()#一行目はseqidとwindow sizeを示すメタ情報なのでスキップ
        for line in f:
            row = line.rstrip().split(",")
            positions.append(int(row[0]))#ポジションは整数に変換 
            values.append(float(row[1]))#GC含量は少数(float)に変換
    fig = plt.figure(figsize=(6,2))
    ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    ax.plot(positions,values)
    fig.savefig(sys.argv[1].replace(".txt",".pdf"),bbox_inches="tight")
````
少し解説。

**with open(sys.argv[1],"r") as f**  
ここで使った```sys.argv[1]```というのは、コマンドライン引数とか呼ばれるやつである。例えば、
````
python plot.py sample_GCratio.txt
````
としてやると、```sys.argv[1]```には```sample2_GCratio.txt```が入る。こんな風に```sys.argv```を使えば、スクリプトを書き変えずともコマンドラインから任意の文字列や数値をスクリプトに渡すことができる。ついでにいうと```sys.argv[0]```にはスクリプトファイル名、つまり```plot.py```が入る。

**row = line.rstrip.split(",")**  
ファイルの中身はカンマ区切りになっていて、1列目がポジション。2列目がGC含量になっている。なので、まず```split(",")```をつかって各行をカンマで分割し、2つの要素からなるListに変換している。

**fig = plt.figure(figsize=(4,2))**  
figureを作成するための土台を作っている。

**ax  = fig.add_axes([0.1,0.1,0.8,0.8])**  
実際に、結果を図示する枠をfig上に生成している。figの縦、横を1.0として捉えたときにどのような枠を設定するかというのが[0.1,0.1,0.8,0.8]の部分。最初の0.1,0.1は起点となるx,yの座標。0.8,0.8は枠の横、縦の長さを示す。

**ax.plot(positions,values)**  
折れ線グラフを表示してくれる。細かい設定は自分でググって勉強しよう。

**fig.savefig(sys.argv[1].replace(".txt",".pdf"),bbox_inches="tight")**  
作成したfigureをpdfに保存している。```bbox_inches="tight"```はおまじない。つけといて損することはあまりないので、つけることをお勧めする。

### 課題1
自分が扱う生物種を対象に、GC_slide2.pyとplot.pyを使ってGC含量の移動プロットをしてみよう。

### 課題2
sample2.fastaを対象にGC skewの移動平均をとるスクリプト"GCskew_slide.py"を実行して、作成されたテキストファイルをplot.pyで描画してみよう。実はsample2.fastaはEscherichia coli str. K-12 substr. MG1655（最も一般的な大腸菌）のゲノムなのだが、図から何か特徴に気づいたかな？
