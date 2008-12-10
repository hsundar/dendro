function writeinCformat(mattr,str)
[a,b,c] = size(mattr);
fid = fopen('cformatfile',str);
for k = 1:c
for i = 1:a
    if (i ==1)
        fprintf(fid,'{\n');
    end
    for j = 1:b
        if(j==1)
            fprintf(fid,'{');
        end
        if(j < b)
            fprintf(fid,'%.15f,',mattr(i,j,k));
        else
            if(i<a)
                fprintf(fid,'%.15f},',mattr(i,j,k));
            else
                fprintf(fid,'%.15f}',mattr(i,j,k));
            end
        end
    end
    fprintf(fid,'\n');
    if(i==a)
        fprintf(fid,'},');
    end
end
fprintf(fid,'\n');
end
fclose(fid);
