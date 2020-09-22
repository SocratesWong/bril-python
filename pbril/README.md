PBril: A Python Based Bril Optimizer
=========================================================



Dead Code Elimiation
-----------------


Execution Example:


		PATH=$PATH:`yarn global bin` bril2json < examples/test/tdce/reassign.bril | python3 pbril/complier.py | bril2txt


Recived Output:


		@main {
		  a: int = const 42;
		  print a;
		}
		
Input:

		cat  examples/test/tdce/reassign.bril
		@main {
		  a: int = const 100;
		  a: int = const 42;
		  print a;
		}
		
As it is seen that the a variable that is being reused is removed as it is dead code.  With other examples, if the variable is not being used in the end of the function DCE will also elimiate it.

LVN Immplemtation
-----------------


LVN is immplemented and it is able to preform all three of the operations.  Should work with the tricker examples.


###Constant Folding


Execution Example:

		PATH=$PATH:`yarn global bin` bril2json < examples/test/tdce/combo.bril | python3 pbril/complier.py | bril2txt
		
Recived Output:

		@main {
		  d: int = const 4;
		  print d;
		}

Input:

		cat examples/test/tdce/combo.bril
		@main {
		  a: int = const 1;
		  b: int = const 2;
		  c: int = add a b;
		  b: int = const 3;
		  d: int = add a b;
		  print d;
		}
		
		
The complier is able to elemiate any constant at compile time.  Instead of having to add the numbers at run time, it replaced it by loading in d directly


###Algebraic Identities (Disable updateinit)


Execution Example:

		PATH=$PATH:`yarn global bin` bril2json < examples/test/lvn/commute.bril | python3 pbril/complier.py | bril2txt
Recived Output:

		@main {
		  a: int = const 4;
		  b: int = const 2;
		  sum1: int = add a b;
		  prod: int = mul sum1 sum1;
		  print prod;
		}

Input:

		cat examples/test/lvn/commute.bril
		# ARGS: -c
		# (a + b) * (b + a)
		@main {
		  a: int = const 4;
		  b: int = const 2;
		  sum1: int = add a b;
		  sum2: int = add b a;
		  prod: int = mul sum1 sum2;
		  print prod;
		}
		
		
The complier understands community laws of addition and mutiplication, and is able to treat tem as the same numbering.  


###Copy propagation


Execution Example:

		PATH=$PATH:`yarn global bin` bril2json < examples/test/lvn/idchain.bril | python3 pbril/complier.py | bril2txt
		
Recived Output:

		@main {
		  x: int = const 4;
		  copy3: int = id x;
		  print copy3;
		}


Input:

		cat examples/test/lvn/idchain.bril
		@main {
		  x: int = const 4;
		  copy1: int = id x;
		  copy2: int = id copy1;
		  copy3: int = id copy2;
		  print copy3;
		}

The complier is able to copy propagation and avoid a needless assignment to copy2.  


