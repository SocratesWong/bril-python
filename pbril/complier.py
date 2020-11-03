import sys
import json
import itertools
from collections import OrderedDict

licmnum=1000;	

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

def dce(fun,livevar):
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
							sys.stderr.write("overwritten deleating {},{}\n".format(db,di));
					lookback[dest]=(b,i);
					unused.add(dest)
					
					
				#if dest in unused
		
		
		for i in unused: 
				db, di = lookback[i]
				if 'label' in blocks[db][0]:
					if i not in livevar[blocks[db][0]['label']]:
						delet[db].add(di);
						done=False;
						changes=changes+1;
						sys.stderr.write(i+" unused deleating {},{}\n".format(db,di));
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
				'args':[self.conical]
			});
		elif instr['op'] =='id':	
			instr.update({
				'op':'id',
				'args':[self.conical]
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
callnum=100;
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
		global callnum
		if op =='call':
			#raise ValueError('No key to gen for Call All keys should be unique for function calls')	
			rt= ('call',callnum);
			callnum=callnum+1;
			return rt;
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
				sys.stderr.write("args is : "+str(instr['args'])+"**********\n");
				instr['args'][0]=carg1; #name update, conical replacement
				if key not in self.lookuptable.keys():
					raise ValueError("We got an BIG PROBLEM ID on undefined value!!!!");
				else:
					
					self.lookuptable[key].addalais(instr['dest'])  
					self.lookuptable[(dest,)]=self.lookuptable[key];
			elif instr['op'] =='call':
				if 'args' in instr:
					instr['args'] = [self.lookupvar(n) for n in arglist];  #conical replacement
				newent= lvn_entery(self,instr['dest'],instr['op'],lastwrite);
				self.lookuptable[newent.conical]=newent;
				key = self.genkey(instr['op']);
				self.lookuptable[key]=newent;
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
					self.lookuptable[key].updatecp(self, instr);
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
		
		
		
def lvn(fun,livevar):
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
	return dce(fun,livevar)
def union(list):
	out=set()
	for i in list:
		out.update(i)
	return out;
def intersect(list):
    if not list:
        return set();
    #sys.stderr.write("Intersecting:"+str(list)+"\n");
    out=set(list[0]);
    for i in list:
        out=set.intersection(out,i)
    return out;
def definedtrans(block,inval):
	gen , kill, used= gku(block);
	return inval.union(gen);
def livetrans(block,inval):
	gen , kill, used= gku(block);
	return used.union(inval-gen);
def gku(block):
	gen =set();
	kill= set()
	used= set();
	for i in block:
		for arg in i.get('args',[]):
			if arg not in gen:
				used.add(arg)
		if 'dest' in i:
			if i['dest'] not in used:
				kill.add(i['dest'])
			gen.add(i['dest'])
		
	return (gen, kill,used);
def succ(instr,next):
	#sys.stderr.write("Succ "+str(instr)+"\n");
	if instr['op'] in ('jmp','br'):
		return instr['labels'];
	if instr['op'] in ('ret',):
		return [];
	if next is None:
		return [];
	return [next];
def appendifexist(dic,name,obj):
	#if not obj:
	#	return
	try:
		dic[name].append(obj)
	except:
		dic[name]=[obj];
class fun_block():
	def __init__(self,fun):
		#sys.stderr.write("Fun are "+ str(fun)+"\n");
		self.bct=0;
		self.dict=OrderedDict();
		wbdict={};
		block=form_blocks(fun["instrs"])
		#sys.stderr.write("block are "+ str(block)+"\n");
		blocks = list(block);
		for x,bk in enumerate(blocks):
			try:
			
				self.dict[bk[0]['label']]=bk[1:]
			except:
				
				self.dict['genlabel.'+str(self.bct)]=bk
				wbdict[x] = 'genlabel.'+str(self.bct)
				#bk.insert(0,{'label':'genlabel.'+str(self.bct)})
				self.bct=self.bct+1;
		for x, val in wbdict.items():
			blocks[x].insert(0,{'label': str(val)})
		fun["instrs"]=flatten(blocks)
		#sys.stderr.write("Dicts are "+ str(dict(self.dict))+"\n");
		#block[:]=fb
	def slover( self,foward,merger, transfer):
		start = list(self.dict.keys())[ 0 if foward else -1 ]  
		if foward:
			in_ed,out_ed= self.edge()
			#sys.stderr.write("in_ed values"+str(in_ed)+"\n");
		else:
			out_ed, in_ed = self.edge()
		inp ={start:set()} 
		out ={bk: set() for bk in self.dict}
		#sys.stderr.write("Keys are "+ str(self.dict.keys())+"\n");
		stack = list(self.dict.keys())
		while stack:
			cur = stack.pop();
			#sys.stderr.write("DEBUG:in_end:"+str(in_ed)+"\nout_ed"+str(out_ed)+"\nout:"+str(out)+"\n");
			#sys.stderr.write("cur:"+str(cur)+"\n");
			inp[cur]=merger(out[p] for p in in_ed[cur])
			out_temp = transfer(self.dict[cur],inp[cur])
			if out[cur] != out_temp:
				stack+=out_ed[cur]
			out[cur]=out_temp
			
		return (inp,out) if foward else (out,inp)
		
	def edge(self):
		pred={name:[] for name in self.dict};
		suc={name:[] for name in self.dict};
		#sys.stderr.write(str(self.dict)+"\n");
		for ct,value in enumerate(self.dict.items()):
			name,block =value;
			try:
				next, _ =list( self.dict.items())[ct+1]
			except:
				sys.stderr.write("Name failed to fetch next block["+name+"] is index "+ str(ct)+"\n")
				#next, _ = self.dict.items()[ct+1]
				next= None;
			#if next:
			#	print("next is "+ next+"\n")
			#else:
			#	print("next is none"+"\n")
			if block:
				for s in succ(block[-1],next):
					#sys.stderr.write("Adding value"+s+"\n");
					appendifexist(suc,name,s);
					appendifexist(pred,s,name);
		return pred,suc;
	def edgelist(self):
		pass
		
	def dom(self):
		
		pred, suc = self.edge();
		dom ={}
		tdom={}
		#sys.stderr.write("Pred list:"+str(pred)+"\n");
		#sys.stderr.write("Suc list:"+str(suc)+"\n");
		if True:
		#for ct,valuet in enumerate(self.dict.items()):
		#	sys.stderr.write("Working on:"+str(valuet[0])+"\n");
			done = False
			dom={}
			for value1 in self.dict.items():
				name,block =value1;
				dom[name]=set([name]);
				dom[name]=set(self.dict.keys());
			#namet,block =valuet;
			#dom[namet]=set(self.dict.keys());
			while not done:
				done = True;
				for value in reversed(self.dict.items()):
					name,block =value;
					newdom= set([name])
					dlist= [];
					for pd in pred[name]:
						dlist.append(dom[pd])
					list = intersect(dlist);
					newdom.update(list);
					done = done & (dom[name] == newdom);
					dom[name]=newdom;
			#tdom[namet]=newdom;

		return dom;
	def strict_dom(self):
		doms=self.dom();
		for name, li in doms.items():
			li.remove(name);
		return doms;
	
	def forntier_dom(self):
		doms=self.dom();
		sdom=self.strict_dom();
		pred,suc = self.edge();
		don_i=map_inv(doms);
		fdoms={}
		ds=set()
		for key in self.dict.keys():
			fdoms[key]=set()
		#for B in self.dict.keys():
		#	for A in self.dict.keys():
		#		pres= set()
		#		pres.update(pred[B]);
		#		prec= pres;
		#		done = False;
				#while not done:
				#	done = True
				#	pres= set()
				#	pres.update(prec);
					#for x in prec:
						#sys.stderr.write("Adding x to prec:"+str(set(pred[x]))+"\n");
						#if not (set(pred[x])):
					#	pres.update(set(pred[x]))
					#if pres!=prec:
						#sys.stderr.write("Comp pres:"+str(pres)+"\n");
						#sys.stderr.write("Comp prec:"+str(prec)+"\n");
					#	done = False
		#			prec=pres
		#		for pre in prec:
		#			if (A not in doms[pre]) and (A in doms[B]) :
		#				fdoms[B].add(A);
		#				sys.stderr.write("A:"+A+" B:"+B+" pre:"+pre+"\n");
		
		#for A in self.dict.keys():
		#	for su in self.dict.keys():
		#		#sucs = suc[B]
		#		B_s = pred[su];
		#		for B in B_s:
		#			if A in doms[B] and A not in doms[su]:
		#				fdoms[B].add(A);
		#				sys.stderr.write("A:"+A+" B:"+B+" su:"+su+"\n");
		for A in doms:
			ds=set()
			for dd in don_i[A]	:
				ds.update(suc[dd])
			for dt in ds:
				if dt not in don_i[A] or dt== A:
					fdoms[A].add(dt)			
		
				
			
		return fdoms
	#def SSA(self):
	#	fontier=self.forntier_dom()
		

def map_inv(inv):
    out = {key: [] for key in inv}
    for tar, intt in inv.items():
        for ink in intt:
            out[ink].append(tar)
    return out
#Imported code from sampsyo/bril.	
def fmt(val):
    """Guess a good way to format a data flow value. (Works for sets and
    dicts, at least.)
    """
    if isinstance(val, set):
        if val:
            return ', '.join(v for v in sorted(val))
        else:
            return '∅'
    elif isinstance(val, dict):
        if val:
            return ', '.join('{}: {}'.format(k, v)
                             for k, v in sorted(val.items()))
        else:
            return '∅'
    else:
        return str(val)
#End of Imported Code.


def to_ssa(fun, livein, liveout,fdom, p,s):
	done= False;
	modified = False;
	blocks = list(form_blocks(fun['instrs']))
	changes =0;
	blockvaldict={}
	valdict={}
	phidict={}
	for block in blocks:
		name=block[0]['label'];
		phidict[name]=[];
		blockvaldict[name]={};
	for bkc, block in enumerate(blocks):
		name=block[0]['label'];
		if not fdom[name]:
			sys.stderr.write("p["+name+"] is "+ str(p[name])+"\n")
			if p[name]:
				#predname= p[name][0];
				blockvaldict[name]=dict(blockvaldict[p[name][0]])
		else:
			sys.stderr.write("creating phi for "+name+"\n");
			for a in livein[name]:
				sys.stderr.write("\t for "+a+"\n");
				varlist=[];
				labellist=[];
			
					
				for ss in p[name]:
					labellist.append(ss);
				
				if a in valdict:
					valdict[a]=a+"."+str(changes)
					blockvaldict[name][a]=a+"."+str(changes);
					changes=changes+1;
				else:
					valdict[a]=a;#+"."+str(changes);
					blockvaldict[name][a]=a;#"."+str(changes);
					#changes=changes+1;
				phiset=(a,varlist,labellist,blockvaldict[name][a])
				phidict[name].append(phiset);
		for i, inst in enumerate(block):
		
			for x,var in enumerate(inst.get('args', [])):
				#update name
				if var in blockvaldict[name]:
					blocks[bkc][i]['args'][x]=blockvaldict[name][var]
				else:
					valdict[var]=var;
					blockvaldict[name][var]=var;
					
				
			if 'dest' in inst:
				a= inst['dest']
				if a in valdict:
					valdict[a]=a+"."+str(changes)
					blockvaldict[name][a]=a+"."+str(changes);
					changes=changes+1;
					blocks[bkc][i]['dest']=blockvaldict[name][a];
				else:
					valdict[a]=a;
					blockvaldict[name][a]=a;
					blocks[bkc][i]['dest']=blockvaldict[name][a];
					
		# update + insert phi
	for bkc, block in enumerate(blocks):
		name=block[0]['label'];
		for phiinsert in phidict[name]:
			for inflow in phiinsert[2]:
				sys.stderr.write("Inflow "+str(inflow)+" Insert "+str(phiinsert[0])+"\n");
				try:
					#sys.stderr.write("blockvaldict "+str(blockvaldict[inflow][phiinsert[0]])+"\n");
					phiinsert[1].append(blockvaldict[inflow][phiinsert[0]])
				except:
					sys.stderr.write("Fallback\n")
					#phiinsert[1].append(inflow)
					pass
			phii={
		        'op': 'phi',
		        'dest': phiinsert[3],
		        'type': 'int',
		        'labels': phiinsert[2],
		        'args': phiinsert[1],
		   	 }
			blocks[bkc].insert(1,phii)
	#sys.stderr.write("**********Function Pass "+str(blocks)+"**********\n");
				
	fun['instrs'] = flatten(blocks)

def licm(fun, livein, liveout,fdom, p,s,fb):
	done= False;
	modified = False;
	blocks = list(form_blocks(fun['instrs']))
	changes =0;
	blockvaldict={}
	valdict={}
	phidict={}
	loopdict={}
	for block in blocks:
		name=block[0]['label'];
		phidict[name]=[];
		#loopdict[name]={};
	for bkc, block in enumerate(blocks):
		name=block[0]['label'];
		loopstart=name;
		loopmembers={name}
		if len(s[name]) !=2:
			continue;
		loopout=s[name][1];
		loopstack=[];
		loopstack.extend(s[name])
		isloop=True;
		while loopstack:
			cname= loopstack.pop();
			if cname is loopout:
				continue
			loopmembers.add(cname);
			if not s[cname]:
				isloop= False;
			
			loopstack.extend(set(s[cname])-loopmembers)
		if isloop:
			#sys.stderr.write("Loop deteched at "+loopstart+" with members"+str(loopmembers)+"exit at "+loopout+"\n");
			loopdict[loopstart]=(bkc,loopstart,loopmembers,loopout);

		else:
			loopstart=name;
			loopmembers={name}
			if len(s[name]) !=2:
				continue;
			loopout=s[name][0];
			loopstack=[];
			loopstack.extend(s[name])
			isloop=True;
			while loopstack:
				cname= loopstack.pop();
				if cname is loopout:
					continue
				loopmembers.add(cname);
				if not s[cname]:
					isloop= False;
			
				loopstack.extend(set(s[cname])-loopmembers)
			if isloop:
				#sys.stderr.write("Loop deteched at "+loopstart+" with members"+str(loopmembers)+"exit at "+loopout+"\n");
				loopdict[loopstart]=(bkc,loopstart,loopmembers,loopout);
				#loopconstant=
	#sys.stderr.write("loopdict "+str(loopdict)+"\n");
	#newinst=[];
	global licmnum
	for t in loopdict:
		if not modified:
			newinst=[];
			
			#sys.stderr.write("Loop deteched at "+str(t)+"\n");
			bn,cloopstart,cloopmem,cloopout= loopdict[t]
			sys.stderr.write("Loop deteched at "+cloopstart+" with members"+str(cloopmem)+"exit at "+cloopout+"\n");
			invarient = livein[cloopstart]	
			killunion= set();
			linst={
			'label': cloopstart
		   	}
							   	
			newinst.append(linst)
			for lm in cloopmem:
				g,k,u= gku(fb.dict[lm]);
				sys.stderr.write("Gen set at  "+lm+" with members"+str(g)+"\n");
				invarient=invarient-g;
				killunion=killunion.union(g);
			sys.stderr.write("Invarient List at   "+lm+" with members"+str(invarient)+"\n");
			#sys.stderr.write("killunion List at   "+lm+" with members"+str(killunion)+"\n");
			
			for block in blocks:
				name=block[0]['label'];
				if name in cloopmem:
					for instr in block:
						#sys.stderr.write("instr   "+str(instr)+"\n");
						if 'op' in instr and instr['op'] in ("add", "mul", 'sub','div','gt','lt','ge','le'):
							arglist=instr.get('args', []);
							if arglist[0] in invarient and arglist[1] in invarient:
								sys.stderr.write("Moving code in   "+lm+" with instr"+str(instr)+"\n");
								linst={
								'op': instr['op'],
								'dest': instr['dest']+'.'+str(licmnum),
								'type': instr['type'],
								'args': list(arglist),
							   	}
							   		
								instr['op']= 'id';
								del instr['args']
								instr['args']= [instr['dest']+'.'+str(licmnum),];
								licmnum=licmnum+1;
								newinst.append(linst);
								modified= True;
			#label rename+ insertblocks
			if modified:
				for block in blocks:
					name=block[0]['label'];
					if name in cloopmem:
						for instr in block:
							if 'labels' in instr:
								instr['labels']= [lbl if  lbl != cloopstart else cloopstart+".old" for lbl in instr['labels']];
							
							if 'label' in instr:
								if name == cloopstart:
									instr['label']=instr['label']+".old";
				for block in blocks:
					name=block[0]['label'];	
					if name == cloopstart+".old":
					#	block[:]= newinst.extend(block)	;	
						for i in reversed(newinst):
							block.insert(0, i)
					
	# update + insert phi
	if False:
		for bkc, block in enumerate(blocks):
			name=block[0]['label'];
			for phiinsert in phidict[name]:
				for inflow in phiinsert[2]:
					sys.stderr.write("Inflow "+str(inflow)+" Insert "+str(phiinsert[0])+"\n");
					try:
						#sys.stderr.write("blockvaldict "+str(blockvaldict[inflow][phiinsert[0]])+"\n");
						phiinsert[1].append(blockvaldict[inflow][phiinsert[0]])
					except:
						sys.stderr.write("Fallback\n")
						#phiinsert[1].append(inflow)
						pass
				phii={
				'op': 'phi',
				'dest': phiinsert[3],
				'type': 'int',
				'labels': phiinsert[2],
				'args': phiinsert[1],
			   	 }
				blocks[bkc].insert(1,phii)
	#sys.stderr.write("**********Function Pass "+str(blocks)+"**********\n");
				
	fun['instrs'] = flatten(blocks)
	return modified;
	
def out_ssa(fun, livein, liveout,fdom, p,s):
	done= False;
	modified = False;
	blocks = list(form_blocks(fun['instrs']))
	changes =0;
	upphi={}
	for block in blocks:
		name=block[0]['label'];
		upphi[name]=[];
		
	for bkc, block in enumerate(blocks):
		name=block[0]['label'];
		for i, inst in enumerate(block):
			if 'op' in inst and inst['op'] is 'phi':
				dest= inst['dest'];
			
				for ct,lbl in enumerate(inst['labels']):
					dat=(dest,inst['args'][ct])
					upphi[lbl].append(dat);
	for bkc, block in enumerate(blocks):
		name=block[0]['label'];
		for phich in upphi[name]:
			iid={
                        'op': 'id',
                        'type': 'int',
                        'args': [phich[1]],
                        'dest': phich[0],
                   	};
			#sys.stderr.write("IID is: "+str(iid)+"**********\n");
			if blocks[bkc][-1]['op'] in TERMINATORS:
				blocks[bkc].insert(-1,iid)
			else: 
				blocks[bkc].insert(len(blocks[bkc]),iid)	
	if True:
		for bkc, block in enumerate(blocks):
			name=block[0]['label'];
			delet=set();
			for i, inst in enumerate(block):
				if 'op' in inst and inst['op'] is 'phi':
					delet.add(i);
		
			block[:]=[inst for i,  inst in enumerate(block) 
				if i not in delet]
	sys.stderr.write("**********Function Pass "+str(blocks)+"**********\n");
	fun['instrs'] = flatten(blocks)

def opt():
	# Code structure adopted from examples, orgrinally located from sampsyo/bril.  
	bril = json.load(sys.stdin)
	change = True;
	cttt=5 
	while change and cttt>0:
		cttt=cttt-1
		change = False;
		if True:
			for fun in bril['functions']:
				done=False;
				passct=0;
				while not done:
					done= True;
					passct=passct+1;
					sys.stderr.write("**********Function Pass "+str(passct)+"**********\n");
					fb=fun_block(fun);
					#sys.stderr.write("Liveness :\n")
					in_,out=fb.slover(False,union,livetrans)
					done= done and dce(fun,out);
					done= done and lvn(fun,out);
					#sys.stderr.write(str(done))
		# DF
		sys.stderr.write("DF :\n")
		for fun in bril['functions']:		
			fb=fun_block(fun);
			#defined
			sys.stderr.write("Defined :\n")
			in_,out=fb.slover(True,union,definedtrans)
			for block in fb.dict:
			    sys.stderr.write('{}:\n'.format(block))
			    sys.stderr.write('  in: '+ fmt(in_[block])+"\n")
			    sys.stderr.write('  out:'+ fmt(out[block])+"\n")
			sys.stderr.write("Liveness :\n")
			in_,out=fb.slover(False,union,livetrans)
			for block in fb.dict:
			    sys.stderr.write('{}:\n'.format(block))
			    sys.stderr.write('  in: '+ fmt(in_[block])+"\n")
			    sys.stderr.write('  out:'+ fmt(out[block])+"\n")
			dom=fb.dom();
			sys.stderr.write("Dom list:"+str(dom)+"\n");
			strict_dom=fb.strict_dom();
			sys.stderr.write("Strict Dom list:"+str(strict_dom)+"\n");
			fdom=fb.forntier_dom();
			sys.stderr.write("Frontier Dom list:"+str(fdom)+"\n");
		if True:
			sys.stderr.write("licm: \n")
			for fun in bril['functions']:
				done=False;
				passct=0;
				if True:
					done= True;
					passct=passct+1;
					sys.stderr.write("**********Function Pass "+str(fun['name'])+"licm**********\n");
					fb=fun_block(fun);
					in_,out=fb.slover(False,union,livetrans);
					fdom=fb.forntier_dom();
					p,s=fb.edge();
					sys.stderr.write("Pre: "+str(p)+"\n");
					sys.stderr.write("Suc: "+str(s)+"\n");
					change=licm(fun,in_, out,fdom,p,s,fb)
					#done= done & ssa(fun)
					#sys.stderr.write(str(done))
					
	if False:
		sys.stderr.write("SSA: \n")
		for fun in bril['functions']:
			done=False;
			passct=0;
			if True:
				done= True;
				passct=passct+1;
				sys.stderr.write("**********Function Pass "+str(fun['name'])+"TO SSA**********\n");
				fb=fun_block(fun);
				in_,out=fb.slover(False,union,livetrans);
				fdom=fb.forntier_dom();
				p,s=fb.edge();
				sys.stderr.write("Pre: "+str(p)+"\n");
				sys.stderr.write("Suc: "+str(s)+"\n");
				to_ssa(fun,in_, out,fdom,p,s)
				#done= done & ssa(fun)
				#sys.stderr.write(str(done))
	if False:
		sys.stderr.write("out SSA: \n")
		for fun in bril['functions']:
			done=False;
			passct=0;
			if True:
				done= True;
				passct=passct+1;
				sys.stderr.write("**********Function Pass "+str("OUT SSA")+"**********\n");
				fb=fun_block(fun);
				in_,out=fb.slover(False,union,livetrans);
				fdom=fb.forntier_dom();
				p,s=fb.edge();
				out_ssa(fun,in_, out,fdom,p,s)
	for fun in bril['functions']:
		done=False;
		passct=0;
		if True:
			pass
	
	json.dump(bril, sys.stdout,indent=2, sort_keys=True)





if __name__ == '__main__':
	opt()
