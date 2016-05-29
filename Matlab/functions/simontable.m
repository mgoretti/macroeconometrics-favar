function  simontable( Matrix, name, label,caption, tablewidth, numbercolumn, namecolumn, headers, note)

% copyright Simon Tièche, mail: tieche.simon@gmail.com 

% LaTeX package required: \usepackage{threeparttable}
% If you have package ctable, remove it. ctable is not of interest. 

% The function simontable( Matrix, name, label,caption, tablewidth, numbercolumn, namecolumn, headers, note)
% creates a simple (no multicolumn) tex table for any number of rows and a
% maximum of 10 columns. Then, to call the table in your main file, you
% must input it. so \input{path/nametable}. If you want to rotate it
% then you can use the package sidewaystable and rotating. Otherwise, do
% not forget that minipage and the swiss-knife tikzipicture exists to put
% the table where you want. 
% 
%
% INPUTS: 
% Matrix - the matrix you want to export in Latex. Note that you may have
% to build that matrix. Loops and concanation may help you.

% name - name/path directory + name that you want to give to your tex
% table. For instance: name = 'Paper/Tables/IC.tex'; puts your table in
% a folder Tables in a folder Paper where Paper is a current folder
% of your current directory and gives the name IC to your tex file.
% 
% label - the label you want to give to your tex table. You have to give one
% if you want to cite the table in your mother document
%
% caption - table title
%
% tablewidth - just give a scalar. For instance, tablewidth = 0.5; will
% create a centered table that fits the middle of the page
%
% numbercolumn - the number and type of columns you want. For instance,
% numbercolum =lcc;
% 
% namecolumn - a row vector of name for your columns name. Warning: each name
% must have the same length --> Hit space button
%
% headers - a scalar that contains the name for the headers and the column
% separation &. For instance, headers = Player & Score 1 & Score 2;
%
% note - the note under the table. Note that if you want to put a ' you
% need to type it twice so ' = '' 

% OUTPUT:
% A beautiful table !



[n m] = size(Matrix);

    
FID = fopen(name, 'w');
fprintf(FID, '\\begin{table}[H] \n \\centering \n');
fprintf(FID, '\\begin{threeparttable} \n');
fprintf(FID, '\\caption{\\textsc{%s}} \\label{%s} \n', caption,label);
fprintf(FID, '\\begin{tabular*}{%s \\textwidth }{%s} \n \\toprule \\toprule  \n', tablewidth, numbercolumn); 
fprintf(FID, ' %s ', headers);
fprintf(FID, '\\tabularnewline \\midrule');
for ii = 1:n
    if m==1
fprintf(FID, ' %s & %8.2f   \\\\ \n', namecolumn(ii,:), Matrix(ii,1));
    elseif m==2
        fprintf(FID, ' %s & %8.2f & %8.2f \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2));
    elseif m ==3
        fprintf(FID, ' %s & %8.2f & %8.2f & %8.2f  \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2), Matrix(ii,3));
    elseif m ==4
        fprintf(FID, ' %s & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2), Matrix(ii,3), Matrix(ii,4));
    elseif m ==5 
        fprintf(FID, ' %s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2), Matrix(ii,3), Matrix(ii,4), Matrix(ii,5));
    elseif m ==6
        fprintf(FID, ' %s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2), Matrix(ii,3), Matrix(ii,4), Matrix(ii,5),  Matrix(ii,6));
    elseif m ==7 
        fprintf(FID, ' %s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2), Matrix(ii,3), Matrix(ii,4), Matrix(ii,5),  Matrix(ii,6), Matrix(ii,7));
    elseif m ==8 
        fprintf(FID, ' %s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2), Matrix(ii,3), Matrix(ii,4), Matrix(ii,5),  Matrix(ii,6), Matrix(ii,7), Matrix(ii,8));
    elseif m == 9
        fprintf(FID, ' %s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2), Matrix(ii,3), Matrix(ii,4), Matrix(ii,5),  Matrix(ii,6), Matrix(ii,7), Matrix(ii,8), Matrix(ii,9));
    elseif m ==10
        fprintf(FID, ' %s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ \n', namecolumn(ii,:), Matrix(ii,1),Matrix(ii,2), Matrix(ii,3), Matrix(ii,4), Matrix(ii,5),  Matrix(ii,6), Matrix(ii,7), Matrix(ii,8), Matrix(ii,9), Matrix(ii,10));
    else 
    end
end
fprintf(FID, ' \\bottomrule \\bottomrule \n \\end{tabular*} \n');
fprintf(FID, '\\begin{tablenotes} \n');
fprintf(FID, '\\small \n');
fprintf(FID, '\\item \\emph{ \\footnotesize{ Notes : %s } } \n',note);
fprintf(FID, '\\end{tablenotes} \n');
fprintf(FID, '\\end{threeparttable} \n');
fprintf(FID, '\\end{table} \n');
fclose(FID);


end



