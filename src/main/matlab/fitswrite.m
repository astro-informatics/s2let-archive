function fitswrite(data, filename, parname, parval)
%Author: R. G. Abraham, Institute of Astronomy, Cambridge University
%        abraham@ast.cam.ac.uk
%Improved by Boris Leistedt bo be compatible with MW and HEALPIX maps

[nrow,ncol]=size(data);
header_cards = [make_card('SIMPLE','T');       ...
      make_card('BITPIX',16);       ...
      make_card('NAXIS',0);          ...
      make_card('EXTEND','T');        ...
      make_card('END');];
  
[nrowbis,ncolbis] = size(header_cards);
n_blanks = 36 - rem(nrowbis,36);
blank_line = setstr(ones(1,80)*32);
header_cards2 = [header_cards; repmat(blank_line,n_blanks,1)];

header_cards2 = [header_cards2;       ...
      make_card('XTENSION','BINTABLE');        ...
      make_card('BITPIX',8);        ...
      make_card('NAXIS',2);          ...
      make_card('NAXIS1',4);      ...
      make_card('NAXIS2',ncol);      ...
      make_card('PCOUNT',0);        ...
      make_card('GCOUNT',1);        ...
      make_card('TFIELDS',1);        ...
      make_card('TFIELDS',1);        ...
      make_card('TTYPE1','SIGNAL');        ...
      make_card('TFORM1','1E');        ...
      make_card('TUNIT1',' ');        ...
      make_card('EXTNAME','BINTABLE');        ...
      make_card(parname, parval);        ...
      make_card('END')];

header_record = make_header_record(header_cards2);
%[ncards,dummy]=size(header_cards);
%fprintf(header_record(1,:));

fid=fopen(filename,'w');
fwrite(fid,header_record','char');

% try to figure out if we need to swap bytes. This is
% imperfect as I don't know the endian-ness of each
% architecture, so I'm only looking for ones I know for 
% sure are big-endian.
friend = computer;
if strmatch(friend,'PCWIN')
   bswap = 'b';
elseif strmatch(friend,'LNX86')
   bswap = 'b';   
elseif strmatch(friend,'ALPHA')
   bswap = 'b';
else
   bswap = 'l';
end

fwrite(fid, data, 'single', 0, 'B');
%fwrite(fid,data,'single',bswap);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function card=make_card(keyword,value)
%MAKE_CARD turns a set of strings into a valid FITS card
%Make keyword field 8 bytes long
lk=length(keyword);
if (lk > 8) & (nargin>1)
	error('Keyword must be less than or equal to 8 characters!')
elseif (lk < 8 )
	keyword=[keyword,setstr(ones(1,8-lk)*32)];
end;

%Deal with both straight keyword and keyword/value pair
if (nargin==1)
	%Keyword without a value
	card=keyword;	
else
	%Key/value pair has an equal sign and space at bytes 9 and 10
	card=[keyword,'= '];

	%Now output the value. The FITS standard wants things to start 
	%in different columns depending on what type of data the
	%value holds, according to the following rules:
	%
	%  Logical: T or F in column 30
	%
	%  Character string: A beginning quote in column 11 and an
	%  ending quote between columns 20 and 80.
	%
	%  Real part of an integer or floating point number: right 
	%  justified, ending in column 30.
	%
	%  Imaginary part: right justified, ending in
	%  column 50, starting after column 30 (NB. I won't bother 
	%  to support an imaginary part in this M-file, and will 
	%  let some radio astronomer who needs it add it if they want).

	if isstr(value)
  	    %Test for logical. If logical it goes in column 30 
		if (length(value)==1) & (strmatch(upper(value),'T') | strmatch(upper(value),'F'))
 			card=[card,setstr(ones(1,19)*32),value];	
		else	
			%Value must be a character string. Pad if less than 8
			%characters long.
			lv=length(value);
		    if (lv > 70)
		   error('Value must be less than 70 characters long!')
		    elseif (lv < 10 )
 	 	   value=[value,setstr(ones(1,8-lv)*32)];
		    end;
			card=[card,'''',value,''''];
		end;	
	else
		%Value must be a number. Convert to a string. Maximum
		%precision is set to 10 digits
		value=num2str(value,10);
		lv=length(value);
	
		%Write this out so it will end on column 30
		card=[card,setstr(ones(1,20-lv)*32),value];	
	end;
end;

%Now pad the output to make it exactly 80 bytes long
card=[card,setstr(ones(1,80-length(card))*32)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hrec=make_header_record(card_matrix)

[nrow,ncol] = size(card_matrix);
n_blanks = 36 - rem(nrow,36);
blank_line = setstr(ones(1,80)*32);
hrec = [card_matrix; repmat(blank_line,n_blanks,1)];
