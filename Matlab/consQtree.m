function nodes = consQtree(maxD)
 figure(1)
 newplot;
hold on
 coor = [0 0;0 1; 1 0; 1 1]
 A = [1 1 1 0;1 1 0 1;1 0 1 1; 0 1 1 1];
 newA = [1 1 0 1 0 0 0 0 0;    
 1 1 1 0 1 0 0 0 0;  
   0 1 1 0 0 1 0 0 0; 
    1 0 0 1 1 0 1 0 0;  
   0 1 0 1 1 1 0 1 0; 
    0 0 1 0 1 1 0 0 1; 
    0 0 0 1 0 0 1 1 0;  
   0 0 0 0 1 0 1 1 1;  
   0 0 0 0 0 1 0 1 1 ]; 
gplot(A,coor,'-k') 
nodes=[0];
but = 1; 
%left-click to create the tree
 while(but == 1)
     [xin, yin, but]=ginput(1)
     %find the leaf that contains this point and divide it.
     num=length(nodes);
     for(i=1:num) 
        [xcor,ycor,lev] = getXYL(nodes(i,1),maxD);
         dimen = getLen(lev); 
        xcor = xcor/(2^maxD);
         ycor = ycor/(2^maxD);  
       if((xin>=xcor) && (xin < (xcor+dimen) )&& (yin>=ycor) && (yin < (ycor+dimen) ))
             if(lev < maxD) 
                child = getChildren(nodes(i,1),maxD);
                 %Replace this node with its first child   
              nodes(i,1) = child(1);   
              %Append the remaining children to the list  
               nodes(num+1,1) = child(2);      
           nodes(num+2,1) = child(3);    
             nodes(num+3,1) = child(4);      
           %Form new Coordinates    
             [newx,newy,dum] = getXYL(child(1),maxD);    
             newC(1,:) = [newx, newy]/(2^maxD);  
                [newx,newy,dum] = getXYL(child(2),maxD);     
            newC(2,:) = [newx, newy]/(2^maxD);   
               [newx,newy,dum] = getXYL(child(3),maxD);     
            newC(4,:) = [newx, newy]/(2^maxD);  
                [newx,newy,dum] = getXYL(child(4),maxD);    
             newC(5,:) = [newx, newy]/(2^maxD);    
              newDim = getLen(lev+1);   
               newC(3,:) = newC(2,:) + [newDim,0];  
               newC(6,:) = newC(5,:) + [newDim,0];   
              newC(7,:) = newC(4,:) + [0,newDim];     
            newC(8,:) = newC(5,:) + [0,newDim];     
            newC(9,:) = newC(8,:) + [newDim,0];      
           gplot(newA,newC,'-k');      
                       end     
        break;    
     end  
   end 
end  

nodes = sort(nodes); 

 %right-click to generate morton ids. 
% while(but == 3)
 %     [xin, yin, but]=ginput(1)
 %  
   %find the leaf that contains this point and print its id. 
%     [dum, num]=size(nodes); 
%     for(i=1:num) 
%         [xcor,ycor,lev] = getXYL(nodes(i),maxD); 
%         dimen = getLen(lev);    
     %         xcor = xcor/(2^maxD);
 %         ycor = ycor/(2^maxD); 
%         if(xin>=xcor && xin <xcor+dimen && yin>=ycor && yin < ycor+dimen ) 
%             text(xcor+(dimen/2),ycor+(dimen/2),num2str(nodes(i))); 
%             break;
 %         end 
%     end 
% end
  hold off;  
