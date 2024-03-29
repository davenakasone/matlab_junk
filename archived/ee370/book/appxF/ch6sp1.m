% Nise, N.S. 
% Control Systems Engineering, 8th ed. 
% John Wiley & Sons, Hoboken, NJ, 07030
%
% Control Systems Engineering Toolbox Version 8.0 
% Copyright � 2019 by John Wiley & Sons, Inc.
%
% Chapter 6: Stability
%
% ch6sp1  (Example 6.2)    MATLAB's Symbolic Math Toolbox may
% be used conveniently to calculate  the values in a Routh table.
% The toolbox is particularly useful for more complicated tables,
% where symbolic objects, such as epsilon, are used. In this example 
% we represent each row of the Routh table by a vector. Expressions are 
% written for subsequent row elements by using the equations given in Table 6.2
% of the text. The MATLAB command det(M) is used to find the determinant 
% of the square matrix, M, as shown for each row element in Table 6.2.
% Further, we test the previous row's first element to see if it is zero.
% If it is zero, it is replaced by epsilon, e, in the next row's calculation.
% The preeceding logic is performed using MATLAB's IF/ELSE/END as shown in the code 
% below. 
% We now demonstrate the making of a Routh table using the Symbolic Math Toolbox 
% or a problem that requires the epsilon method to complete the table. The following 
% program produces the Routh table for Example 6.2 in the text. Also, for clarity, 
% we convert all rows to symbolic objects, simplify, and pretty print after forming 
% the table. CAUTION: In general, the results of this program are not valid if an 
% entire row is zero as e approaches zero, such as [e 0 0 0]. This case must be 
% handled differently, as discussed in text Section 6.3 in the subsection, 
% "Entire Row is Zero."

'(ch6sp1)  Example 6.2'       % Display label.
% -det([si() si();sj() sj()])/sj()
                              % Template for use in each cell.
syms e                        % Construct a symbolic object for 
                              % epsilon.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s5=[1 3 5 0 0];               % Create s^5 row of Routh table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s4=[2 6 3 0 0];               % Create s^4 row of Routh table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if -det([s5(1) s5(2);s4(1) s4(2)])/s4(1)==0
	s3=[e...
 -det([s5(1) s5(3);s4(1) s4(3)])/s4(1)   0   0];
                              % Create s^3 row of Routh table 
                              % if 1st element is 0.
else
    s3=[-det([s5(1) s5(2);s4(1) s4(2)])/s4(1)...
 -det([s5(1) s5(3);s4(1) s4(3)])/s4(1)   0   0];
                              % Create s^3 row of Routh table 
                              % if 1st element is not zero.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if -det([s4(1) s4(2);s3(1) s3(2)])/s3(1)==0
	s2=[e ... 
 -det([s4(1) s4(3);s3(1) s3(3)])/s3(1)   0   0];
                              % Create s^2 row of Routh table 
                              % if 1st element is 0.
else	
    s2=[-det([s4(1) s4(2);s3(1) s3(2)])/s3(1) ... 
 -det([s4(1) s4(3);s3(1) s3(3)])/s3(1)   0   0];
                              % Create s^2 row of Routh table 
                              % if 1st element is not zero.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if -det([s3(1) s3(2);s2(1) s2(2)])/s2(1)==0
	s1=[e ... 
 -det([s3(1) s3(3);s2(1) s2(3)])/s2(1)   0   0];
                              % Create s^1 row of Routh table 
                              % if 1st element is 0.
else
s1=[-det([s3(1) s3(2);s2(1) s2(2)])/s2(1) ... 
 -det([s3(1) s3(3);s2(1) s2(3)])/s2(1)   0   0];
                              % Create s^1 row of Routh table 
                              % if 1st element is not zero
end							  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%							  
s0=[-det([s2(1) s2(2);s1(1) s1(2)])/s1(1) ... 
 -det([s2(1) s2(3);s1(1) s1(3)])/s1(1)   0   0];
                              % Create s^0 row of Routh table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
's5'                          % Display label. 
s5=sym(s5);                   % Convert s5 to a symbolic object.
s5=simplify(s5);              % Simplify terms in s^5 row.
pretty(s5)                    % Pretty print s^5 row.
's4'                          % Display label. 
s4=sym(s4);                   % Convert s4 to a symbolic object.
s4=simplify(s4);              % Simplify terms in s^4 row.
pretty(s4)                    % Pretty print s^4 row.
's3'                          % Display label. 
s3=sym(s3);                   % Convert s3 to a symbolic object.
s3=simplify(s3);              % Simplify terms in s^3 row.
pretty(s3)                    % Pretty print s^3 row.
's2'                          % Display label.
s2=sym(s2);                   % Convert s2 to a symbolic object.
s2=simplify(s2);              % Simplify terms in s^2 row.
pretty(s2)                    % Pretty print s^2 row.
's1'                          % Display label.
s1=sym(s1);                   % Convert s1 to a symbolic object.
s1=simplify(s1);              % Simplify terms in s^1 row.
pretty(s1)                    % Pretty print s^1 row.
's0'                          % Display label.
s0=sym(s0);                   % Convert s0 to a symbolic object.
s0=simplify(s0);              % Simplify terms in s^0 row.
pretty(s0)                    % Pretty print s^0 row.
