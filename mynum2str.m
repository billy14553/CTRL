function str = mynum2str(num)
bit = 4;
num1 = num;
while num1>=10
    num1 = num1/10;
    bit = bit -1;
end
strMath = "%."+num2str(bit)+"f";
str = sprintf(strMath,num);
end

