## DNA sequence Statistics wiith Python

ここでは、BioPythonを使ってゲノム配列データを解析する手法について扱う。BioPythonはめちゃくちゃ機能が多いんだけど、マジで解説が少ない。基本機能の紹介しつつ、Biopythonを使った場合の書き方と使わない場合の書き方を書いておく。 以下にサンプルファイルが用意してあるので、自分の作業ディレクトリの直下にコピーすること。

#### Fastaの読み込み->反復配列->GC*ratio->GC*skew

**with Biopython**
````Python 
import os 
import sys 
from Bio import SeqIO 
from Bio import SeqUtils 
if __name__ == "__main__": 
	for fna in SeqIO.parse("sample.fasta","fasta") 
		print(fna.id)
		print(fna.seq) 
		print(fna.seq.reverse_complement()) 
		print(SeqUtils.GC(fna.seq)) 
		print(SeqUtils.GC_skew(fna.seq))`
````

**without Biopython**

\````

\````