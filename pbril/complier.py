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
		delet=[];
		#sys.stderr.write("JSON dump\n")
		#sys.stderr.write(str(blocks));
		for b, block in enumerate(blocks):
			delet.append(set())
			for i, inst in enumerate(block):
				#sys.stderr.write("inside string\n")
				#sys.stderr.write(str(inst))
				#print(unused);
				for var in inst.get('args', []):
					#print(var);
					removeifexist(unused,var);
				
				if 'dest' in inst:
					dest=inst['dest'];
					if dest in unused:
						db, di = lookback[dest]
						if db == b:
							delet[db].add(di);
							#sys.stderr.write(str(type(lookback[dest])));
							#sys.stderr.write(str(lookback[dest]));
							#sys.stderr.write("\n");
							done=False;
							changes=changes+1;
							#sys.stderr.write("overwritten deleating {},{}\n".format(db,di));
					lookback[dest]=(b,i);
					unused.add(dest)
					
					
				#if dest in unused
		
		
		for i in unused:
			db, di = lookback[i]
			delet[db].add(di);
			done=False;
			changes=changes+1;
			#sys.stderr.write(i+" unused deleating {},{}\n".format(db,di));
		for b,block in enumerate(blocks):
			block[:]=[inst for i,  inst in enumerate(block) 
				if i not in delet[b]]
		#sys.stderr.write("Iteration done\n");
		
	fun['instrs'] = flatten(blocks)
	sys.stderr.write("DCE pass done optimization changes: "+str(changes)+"\n");
	return changes == 0;

def canon(value):
	return  Value(value.op, tuple(sorted(value.args))) if value.op in ('add', 'mul') else value;

# Adopted from examples, orgrinally located from sampsyo/bril.  
FOLDABLE_OPS = {
    'add': lambda a, b: a + b,
    'mul': lambda a, b: a * b,
    'sub': lambda a, b: a - b,
    'div': lambda a, b: a // b,
    'gt': lambda a, b: a > b,
    'lt': lambda a, b: a < b,
    'ge': lambda a, b: a >= b,
    'le': lambda a, b: a <= b,
} 
# end of Imported Code
class lvn_entery:
	
	def __init__(self,manger,op,var): # input entery
		self.num=manger.lvnnum;
		self.input=True;
		self.const=False;
		self.exp= False;
		self.call= False;
		manger.lvnnum=manger.lvnnum+1;
		self.conical=var;
		self.alais=[var];
		self.eq=(var,);
		
		
	def __init__(self,manger,outvar,op='null',lastwrite=None, constval=None, arg1=None, arg2=None):
		self.num=manger.lvnnum;
		self.input=False;
		self.const=False;
		self.exp= False;
		self.call= False;
		if op =='null':
			manger.lvnnum=manger.lvnnum+1;
			var=outvar;
			self.conical=var;
			self.alais=[var];
			self.eq=(var,);
		else :
			if op =='call':
				self.call=True;
				self.eq=('call',self.num);
			elif op =='const':
				self.const=True;
				self.constval=constval;
				self.eq=('const',self.constval);
			elif op =='id':
				raise ValueError('Id should not generate an Entry.  It should be added to existing entry!!!')	
			elif op in ("add", "mul", 'sub','div','gt','lt','ge','le'):
				self.exp=True;
				self.eq=(op,arg1,arg2);
				self.arg1=arg1;
				self.arg2=arg2;
			else:
				raise ValueError("Undefined Operation!!!!");
			manger.lvnnum=manger.lvnnum+1;
			if lastwrite:
				self.conical=outvar;
			else:
				self.conical="lvn.{}".format(self.num);
			self.alais=[outvar];
		
	def addalais(self, var):
		self.alais.append(var);
	def updatecp(self, manger, instr):
		if instr['op'] =='call': 
			raise ValueError('Should not happen');
		elif instr['op'] =='const':
			#del instr['args']
			instr.update({
				'op':'id',
				'args':self.conical
			});
		elif instr['op'] =='id':	
			instr.update({
				'op':'id',
				'args':self.conical
			});
		elif instr['op'] in ("add", "mul", 'sub','div','gt','lt','ge','le'):
			instr.update({
				'op':'id',
				'args':[self.conical],
			});
		
	def updateinit(self, manger, instr):
		#return
		#val=tmp;
		if instr['op'] in ("add", "mul", 'sub','div','gt','lt','ge','le'):
			#sys.stderr.write("KEYS dump"+str(manger.lookuptable.keys())+"\n");
			#sys.stderr.write("KEYSA dump"+str((self.arg1,))+"\n");
			#sys.stderr.write("Folding check: "+str((self.arg1,) in manger.lookuptable.keys())+","+str((self.arg2,) in manger.lookuptable.keys())+","+"\n");
			if (self.arg1,) in manger.lookuptable.keys() and (self.arg2,) in manger.lookuptable.keys() and manger.lookuptable[(self.arg1,)].const and manger.lookuptable[(self.arg2,)].const:
				try:
					self.constval=FOLDABLE_OPS[instr['op']](manger.lookuptable[(self.arg1,)].constval,manger.lookuptable[(self.arg2,)].constval);
					self.const=True;
					self.exp= False;
					del instr['args']
					instr.update({
						'op':'const',
						'value':self.constval,
					});
				except:
					pass;
		else:
			pass
class Lvn_scope_manager:
	
	
	
	def __init__(self):
		self.num=0;
		self.lookuptable={}
		self.lvnnum=0;
		self.lvntb={};
		
	def lookupvar(self, var):
		if (var,) not in self.lookuptable.keys():
			newent= lvn_entery(self,var,'null');
			self.lookuptable[(var,)]=newent;
			self.lookuptable[newent.conical]=newent;
			return newent.conical;
		else:
			return self.lookuptable[(var,)].conical;
			
	def genkey(self, op, carg1=None,carg2=None,value=None):
		if op =='call':
			raise ValueError('No key to gen for Call All keys should be unique for function calls')	
		elif op =='const':
			return ('const',carg1);
		elif op =='id':
			return (carg1,);
			#raise ValueError('No key to gen for ID calls')	
		elif op  in ("add", "mul"):
		
			return (op,tuple(sorted([carg1,carg2])));
		elif op in( 'sub','div','gt','lt','ge','le'):
			return (op,(carg1,carg2));
		else:
			#print(op)
			raise ValueError("Undefined OP"+str(op)+" in Key Operation!!!!");
	
	def processinst(self, instr, lastwrite):
		#changes=0;
		arglist=instr.get('args', []);
		if 'dest' not in instr: 
			#nothing to be done
			pass
		else:
			dest=instr['dest'];
			if instr['op'] =='id': 
				carg1=self.lookupvar(arglist[0])
				key = self.genkey(instr['op'],carg1);
				#lvnentery=genkey(instr['op'],carg1)
				instr['args'][0]=carg1; #name update, conical replacement
				if key not in self.lookuptable.keys():
					raise ValueError("We got an BIG PROBLEM ID on undefined value!!!!");
				else:
					
					self.lookuptable[key].addalais(instr['dest'])  
					self.lookuptable[(dest,)]=self.lookuptable[key];
			elif instr['op'] =='call':
				if 'args' in instr:
					instr['args'] = [lookupvar(n) for n in arglist];  #conical replacement
				newent= lvn_entery(self,instr['dest'],instr['op'],lastwrite);
				self.lookuptable[newent.conical]=newent;
				self.lookuptable[key].addalais(dest)  #conical replacement
				self.lookuptable[(dest,)]=self.lookuptable[key];
			elif instr['op'] =='const':
				carg1=instr['value'];
				# nothing to update for arguments
				key = self.genkey(instr['op'],carg1);
				if key not in self.lookuptable.keys():
					newent= lvn_entery(self,instr['dest'],instr['op'],lastwrite,carg1);
					self.lookuptable[newent.conical]=newent;
					self.lookuptable[key]=newent;
					
				else:
					lookuptable[key].updatecp(self, instr);
				self.lookuptable[key].addalais(dest);
				self.lookuptable[(dest,)]=self.lookuptable[key];
			elif instr['op']  in ("add", "mul", 'sub','div','gt','lt','ge','le'):
				carg1=self.lookupvar(arglist[0])
				carg2=self.lookupvar(arglist[1])
				instr['args'][0]=carg1; #name update, conical replacement
				instr['args'][1]=carg2; #name update, conical replacement
				key = self.genkey(instr['op'],carg1,carg2);
				if key not in self.lookuptable.keys():
					newent= lvn_entery(self,instr['dest'],instr['op'],lastwrite,None,carg1,carg2);
					self.lookuptable[newent.conical]=newent;
					self.lookuptable[key]=newent;
					sys.stderr.write("Atempting Folding: "+str(instr['op'])+","+str(carg1)+","+str(carg2)+"\n");
					newent.updateinit(self,instr);
				else:
					self.lookuptable[key].updatecp(self,instr);
				self.lookuptable[key].addalais(dest);
				self.lookuptable[(dest,)]=self.lookuptable[key];
			
#def isvarused(blocks, cindex,var):
#	for n in range(cindex,len(blocks)):
		
		
		
def lvn(fun):
	done= False;
	modified = False;
	blocks = list(form_blocks(fun['instrs']))
	changes =0;
	while not done:
		done=True;
		#used=set();
		for block in blocks:
			num=0;
			var2num={}; #value tbl
			#conical={}
			num2const={};
			num2val={}
			outlines=[False]*len(block);
			#valswap={}; changed to val swap class
			seen=set();
			delet=set();
			manger=Lvn_scope_manager();
			for i, inst in reversed(list(enumerate(block))):
				if 'dest' in inst: 
					dest=inst['dest'];
					outlines[i]= True if dest not in seen else False;
			for i, inst in enumerate(block):
				manger.processinst(inst,outlines[i]);
			
				
	fun['instrs'] = flatten(blocks)
	sys.stderr.write("LVN pass done Running DCE to clean up: \n");
	#return changes == 0;
	return dce(fun)

#Imported code from sampsyo/bril.
def opt():
	# Code structure adopted from examples, orgrinally located from sampsyo/bril.  
	bril = json.load(sys.stdin)
	
	for fun in bril['functions']:
		done=False;
		passct=0;
		while not done:
			done= True;
			passct=passct+1;
			sys.stderr.write("**********Function Pass "+str(passct)+"**********\n");
			done= done and dce(fun);
			done= done and lvn(fun);
			#sys.stderr.write(str(done))
			
	json.dump(bril, sys.stdout,indent=2, sort_keys=True)




#End of Imported Code.
if __name__ == '__main__':
	opt()