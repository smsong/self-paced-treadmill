function [ bytes ] = int16toBytes( numbers )
%int16toBytes Gets a vector (Mx1) of numbers in the -2^15:2^15-1 interval,
%and returns an array of twice the size(Mx2), where every element is one byte  of
%the 2 byte representation of the original number (taken as a int16). Bytes
%are returned as an integer in the 0-255 range
%AUTHOR:
%Pablo Iturralde - pai7@pitt.edu
%Last update: Feb 21st 2013 - 12:00

N=length(numbers);
bytes=zeros(N,2,'uint8');
for i=1:N
%Check if numbers are within range, if not, throw a warning
    if numbers(i)>(2^15-1)
        disp(['Warning: number out of range for conversion. Saturating entry #' num2str(i)])
        numbers(i)=2^15-1;
    elseif numbers(i)<(-2^15)
        disp(['Warning: number out of range for conversion. Saturating entry #' num2str(i)])
        numbers(i)=-2^15;
    end
    if numbers(i)<0 %Negative numbers
        aux=2^15+numbers(i); %aux is the 15-bit 2's complement of the number
        byte1=128; %Sign is stored in byte1 (MSB)
    else %Positive numbers
        aux=numbers(i); %aux is just the number
        byte1=0; %Sign is stored in byte1 (lsb)
    end
        byte1=byte1+uint8(floor(aux/2^8)); %First byte is the sign (already stored) plus the conversion of the MSB
        byte2=uint8(aux-2^8 * floor(aux/2^8)); %Conversion of the lSB, by substraction of the MSB

    bytes(i,1)=byte1;
    bytes(i,2)=byte2;
end


end

