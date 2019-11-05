function sendTreadmillPacket(payload,t)


fwrite(t,payload,'uint8');
%Read reply (no need)
% read_format=fread(t,1,'uint8');
% speeds=fread(t,4,'int16');
% incl=fread(t,1,'int16');
% padding=fread(t,21,'uint8');
end

