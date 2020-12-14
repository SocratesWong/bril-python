This is a shell base immplemnation of a trace-based [very naive] JIT complier.  
Usage is jit.sh [inputfile].  
Outputs in out.txt in bril JSON format.  It is undefined behivor if a function is being invoked in it's execution path.
Results:
examples/test/df/fact.bril
Regular dynamic instructions: 62
Jit dynamic instuctions: 58
benchmarks/loopfact.bril [Trace Args= 30, input args =30]:
Regular dynamic instructions: 402
Jit dynamic instuctions: 376
[Trace Args= 30, input args =20]:
Regular dynamic instructions: 272
Jit dynamic instuctions: 522

The JIT complier works, but is it very naive and the recover protocal is relative simple.  It does provide a preformance benfit if the path taken is the same as the trace, and will always return the correct results in the absance of undefined behavior. 
