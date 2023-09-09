%-----Initialization-----%
r = 8;
n = 15;
k = 11;
dlen = 22;
ch_snr = 30;
snr_decrease = 1;
alphabet = 1:5;
prob = [0.25 , 0.25 , 0.2 , 0.2 , 0.1];
sourc_code_algorithm = input("choose your source coding algorithm" +newline+
"arithmetic(1) huffman(2)"+newline);
loop_src = 100;
loop_ch_snr = 1;
ser_vector = zeros(1, loop_ch_snr*loop_src);
ber_vector = zeros(1, loop_ch_snr*loop_src);
ber_mean_snr = zeros(1,loop_ch_snr);
ser_mean_snr = zeros(1,loop_ch_snr);

for i = 1:loop_ch_snr
	for j = 1:loop_src
		%-----Source-----%
		source = randsrc(1, dlen,[alphabet; prob]);
		seq = source;
		counts = 100*prob;
		%-----source encoding-----%
		if sourc_code_algorithm ==1
			%{arithmetic encoder}%
			code = arithenco(seq, counts);
		elseif sourc_code_algorithm ==2
			%{huffman encoder}%
			[dict, avglen] = huffmandict(alphabet, prob);
			code = huffmanenco(source, dict);
		end
		%-----Channel encoder-----%
		ch_enc = encode(code,n,k,'hamming/binary');
		%-----modulation-----%
		ch_enc = ch_enc';
		ch_enc_int = bit2int(ch_enc,log2(r));
		modulated_code = pammod(ch_enc_int,r,0,'Gray');
		%-----Phisycal channel error-----%
		filter_impulse = rcosdesign(0.25,6,4);
		filter_code = upfirdn(modulated_code, filter_impulse);
		noisy_code = awgn(filter_code, ch_snr, 'measured');
		%-----Demodulation-----%
		demodulated_code = pamdemod(real(noisy_code), r, 0, 'Gray');
		demodulated_code_bin = int2bit(demodulated_code, log2(r));
		demodulated_code_bin = demodulated_code_bin';
		%-----Channel decode-----%
		ch_dec = decode(demodulated_code_bin, n, k, 'hamming/binary');
		%-----Source decoding-----%
		if sourc_code_algorithm ==1
			%{arithmetic decode}%
			src_dec = arithdeco(ch_dec, counts, length(seq));
		elseif sourc_code_algorithm ==2
			%{huffman decoder}%
			src_dec = huffmandeco(ch_dec, dict);
		end
		%-----Error rate-----%
		src_dec = src_dec(:, 1:min(length(src_dec), length(source)) );
		source = source(:, 1:min(length(src_dec), length(source)) );
		[number1, ser] = symerr(src_dec, source);
		[number2, ber] = biterr(src_dec, source);
		ber_vector(1,j+( (i-1)*loop_src) ) = ber;
		ser_vector(1,j+( (i-1)*loop_src) ) = ser;
	end
	ch_snr = ch_snr - snr_decrease;
end

for i=1:loop_ch_snr
	ber_mean_snr(i) = mean(ber_vector(1, (i-1)*loop_src +1: i*loop_src ));
	ser_mean_snr(i) = mean(ser_vector(1, (i-1)*loop_src +1: i*loop_src ));
end
figure(3)
subplot(1,2,1)
plot(ber_mean_snr)
set(gca, 'YScale', 'log')
xlabel("Snr(db)")
ylabel("BER(log)")
grid
subplot(1,2,2)
plot(ser_mean_snr)
set(gca, 'YScale', 'log')
xlabel("Snr(db)")
ylabel("SER(log)")
grid