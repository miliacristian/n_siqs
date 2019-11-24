Istruzioni per l'uso:
	Istruzioni per eseguire il programma con debug:
		1)scegliere 2 numeri primi p,q moltiplicarli e trovare n=p*q (es p=17 q=2017 n=17*2017=34289) n.b non esagerare con le cifre di n 				(massimo 80 cifre)
		2)scrivere nel file solamente "n" trovato(niente spazi e niente caratteri oltre al numero)
		3)andare nella directory dove sono presenti i file sorgente e qui aprire un terminale
		4)se non presente il file eseguibile "runtest" digitare nel terminale il comando "gcc -o runtest runtest.c" (senza apici!)
		5)per eseguire digitare nel terminale il comando "./runtest"  (senza apici!)
		6)n.b. :il runtest è molto più lento del run

	Istruzioni per eseguire il programma senza debug:
		1)scegliere 2 numeri primi p,q moltiplicarli e trovare n=p*q (es p=17 q=2017 n=17*2017=34289)
		2)scrivere nel file solamente "n" trovato(niente spazi niente caratteri oltre al numero)
		3)andare nella directory dove sono presenti i file sorgente e qui aprire un terminale
		4)se non presente il file eseguibile "run" digitare nel terminale il comando "gcc -o run run.c" (senza apici!)
		5)per eseguire digitare nel terminale il comando "./run"  (senza apici!)

 


Istruzioni per nerd:
		stringa per compilare il programma:"make install"  (senza apici!)
		stringa per eseguire il programma senza debug valgrind: "./criv_quad number.txt"  (senza apici!)
		stringa per eseguire il programma con debug valgrind: "valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./criv_quad number.txt"       (senza apici!)
