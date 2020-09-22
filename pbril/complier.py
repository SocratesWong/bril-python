import sys
import json
import itertools

#Imported code from sampsyo/bril.
TERMINATORS = 'br', 'jmp', 'ret'
def flatten(ll):
    # Code adopted from examples, orgrinally located from sampsyo/bril.  
    """Flatten an iterable of iterable to a single list.
    """
    return list(itertools.chain(*ll))
    
def form_blocks(instrs):
    # Code adopted from examples, orgrinally located from sampsyo/bril.  
    """Given a list of Bril instructions, generate a sequence of
    instruction lists representing the basic blocks in the program.
    Every instruction in `instr` will show up in exactly one block. Jump
    and branch instructions may only appear at the end of a block, and
    control can transfer only to the top of a basic block---so labels
    can only appear at the *start* of a basic block. Basic blocks may
    not be empty.
    """

    # Start with an empty block.
    cur_block = []

    for instr in instrs:
        if 'op' in instr:  # It's an instruction.
            # Add the instruction to the currently-being-formed block.
            cur_block.append(instr)

            # If this is a terminator (branching instruction), it's the
            # last instruction in the block. Finish this block and
            # start a new one.
            if instr['op'] in TERMINATORS:
                yield cur_block
                cur_block = []

        else:  # It's a label.
            # End the block here (if it contains anything).
            if cur_block:
                yield cur_block

            # Start a new block with the label.
            cur_block = [instr]

    # Produce the final block, if any.
    if cur_block:
        yield cur_block

#End of Imported Code

def removeifexist(s,value):
	if value in s:
		s.remove(value);
def dcepscan(fun, used,overwritten):
	pass
def dce(fun):
	done= False;
	modified = False;
	blocks = list(form_blocks(fun['instrs']))
	changes =0;
	while not done:
		done=True;
		#used=set();
		unused=set();
		lookback={};
		delet=set();
		for block in blocks:
			for i, inst in enumerate(block):
				for var in inst.get('args', []):
					removeifexist(unused,var);
				
				if 'dest' in inst:
					dest=inst['dest'];
					if dest in unused:
						delet.add(lookback[dest]);
						done=False;
						changes=changes+1;
					lookback[dest]=i;
					unused.add(dest)
					
					
				#if dest in unused
		
		
		for i in unused:
			delet.add(lookback[i]);
			done=False;
			changes=changes+1;
		block[:]=[inst for i,  inst in enumerate(block) 
			if i not in delet]
		
	fun['instrs'] = flatten(blocks)
	sys.stderr.write("DCE pass done optimization changes: "+str(changes)+"\n");
	return changes == 0;



#Imported code from sampsyo/bril.
def opt():
	# Code adopted from examples, orgrinally located from sampsyo/bril.  
	bril = json.load(sys.stdin)
	
	for fun in bril['functions']:
		done=False;
		passct=0;
		while not done:
			done= True;
			passct=passct+1;
			sys.stderr.write("**********Function Pass "+str(passct)+"**********\n");
			done= done and dce(fun);
			#sys.stderr.write(str(done))
			
	json.dump(bril, sys.stdout,indent=2, sort_keys=True)




#End of Imported Code.
if __name__ == '__main__':
	opt()
