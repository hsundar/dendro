function nodes = dispOct(nodes,maxD,figNum) 
figure(figNum) 
newplot;
hold on 
 numNodes = length(nodes);
 A = [1 1 1 0; 1 1 0 1;
 1 0 1 1; 0 1 1 1]; 
for(i = 1:numNodes) 	
[xcor,ycor,lev] = getXYL(nodes(i),maxD);
 	dimen = getLen(lev);
 	xcor = xcor/(2^maxD); 
	ycor = ycor/(2^maxD); 	
coor=[xcor, ycor;xcor, (ycor+dimen); 
(xcor+dimen), ycor; 
(xcor+dimen),(ycor+dimen)];	
 	gplot(A,coor,'-k') 
end
 hold off  
