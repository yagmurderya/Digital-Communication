%% teorik
x = 0:1:15;          % x ekseni vektörel olarak 0-15 dB arası tanımlandı
N_0 = (10.^-(x/10)); % 10log_10(x) = 1/N_0'dan N_0'ın değerleri hesaplandı

% a) eşit olasılıklı bitler _______________________________________________
% işlem kolaylığı için Eb=1 seçildi (Eb = (3/4)(A^2)T)
a2t_a = 4/3;

% a_i değerleri a için tanımlandı
a1_a = (3/2)*a2t_a; % a1
a2_a = -a2t_a;      % a2

sgm_02_a = (5/4).*N_0*a2t_a;        % (sigma_0)^2
gama_a = (1/4)*a2t_a;               % gama_0
Ed_a = (5/2)*a2t_a;                 % Ed

Pb_a = qfunc(sqrt(Ed_a./(2.*N_0))); % Pb değeri teorik olarak tanımlandı

% b) P(1) = 1/5 ve P(0) = 4/5 _____________________________________________
% işlem kolaylığı için Eb=1 seçildi (Eb = (3/5)(A^2)T)
a2t_b = 5/3;

% a_i değerleri b için tanımlandı
a1_b = (3/2)*a2t_b; 
a2_b = -a2t_b;

sgm_02_b = (5/4).*N_0*a2t_b;            % (sigma_0)^2
gama_b = (1/2).*N_0*log(4)+(1/4)*a2t_b; % gama_0
Ed_b = (5/2)*a2t_b;                     % Ed

Pb_b = qfunc((a1_b - gama_b)./sqrt(sgm_02_b))*(1/5) + qfunc((gama_b - a2_b)./sqrt(sgm_02_b))*(4/5);

%% simulasyon
N = 10*10^6;                % her SNR değeri için 10 milyon bit üretilecek

% a) eşit olasılıklı bitler________________________________________________
s_a = rand(N, length(N_0)); % s tanımlandı, length(N_0): SNR sayısı
s_a(s_a >= 0.5) = a1_a;     % 1 bitleri için a1
s_a(s_a < 0.5) = a2_a;      % 0 bitleri için a2

noise_a = sqrt(sgm_02_a).*randn(N, length(N_0)); % gürültü tanımlandı
% randn fonksiyonu ile normal dağılıma sahip N x length(N_0) uzunluğunda
% bir değişken tanımlanarak gerekli sigma_0 değeri ile çarpıldı 
% böylelikle varyans = 1*sigma_0 oldu

z_a = s_a + noise_a;        % örneklenmiş z işareti tanımlandı
% s_a ve noisa_a N x length(N_0) uzunluğunda olduğu için z_a vektörü de bu
% uzunluktadır, yani her SNR değeri için 10 milyon z işareti

% karar aşaması
z_a(z_a > gama_a) = a1_a;   % eğer z işareti gama değerinden büyükse a1'e
z_a(z_a < gama_a) = a2_a;   % küçükse a2'ye eşit olacak

err_a = sum(z_a ~= s_a)./N; % karar verildikten sonra elde edilen işaret 
% ile giriş işareti karşılaştırılarak sum() fonksiyonu ile eşit olmadıkları
% noktalarda err_a değeri 1 artıyor. Her SNR değeri için bu işlem yapıldığı
% için, err 1 x length(N_0) uzunluğunda. Daha sonra bit sayısına (N)
% bölünerek hata oranı bulunuyor

% b) P(1) = 1/5 ve P(0) = 4/5 _____________________________________________
s_b = rand(N, length(N_0));   % s tanımlandı
noise_b = sqrt(sgm_02_b).*randn(N, length(N_0)); % a şıkkındakiyle tamamen
% aynı şekilde gürültü oluşturuldu fakat bu şıkta a2t değeri değiştiği için
% sigma_0 değeri de değişeceğinden gürültünün varyansı değişir

s_b_ = zeros(N, length(N_0)); % farklı olasılıklı bitlere göre a1 ve a2 
% ile doldurmak için aynı boyutlu sıfırlardan oluşan yeni s tanımlandı 
% çünkü başlangıçta tanımlanan s_b matrisi kullanıldığında grafikte sapma
% meydana geldi

s_b_(s_b <1/5) = a1_b;        % random üretilen sayıların 1/5den küçük olma 
% olasılığı (rand fonksiyonu düzgün dağılımda 0-1 arası sayı ürettiği
% için) 1/5'e yani P(1)'e eşittir. Öyleyse s = a1 olmalı
s_b_(s_b > 1/5) = a2_b;       % 1/5'den büyük olma olasılığı 4/5 = P(0)'dır 
% yani s = a2 olmalı

z_b = s_b_ + noise_b;         % örneklenmiş z işareti tanımlandı

% karar aşaması
z_b(z_b > gama_b) = a1_b;     % eğer z işareti gama değerinden büyükse a1'e
z_b(z_b < gama_b) = a2_b;     % küçükse a2'ye eşit olacak

err_b = sum(z_b ~= s_b_)./N;  % karar aşamasından sonra farklı olan bitler
% a şıkkındaki gibi hesaplandı. Hata oranını bulmak için bit sayısına 
% bölündü

%% BHO eğrilerinin çizdirilmesi

semilogy(x, Pb_a, 'y', 'LineWidth',1)    % teorik olarak hesaplanan a şıkkı
% çizdirildi
hold on                   % grafikleri üst üste çizdirebilmek için
title('ELM365 - MATLAB Projesi - 171024011') % başlık eklendi
ylim([10^(-6), 10^(0)])   % y ekseni sınırlandırıldı
xlabel('SNR = 10log_{10}(Eb/N0)'), ylabel('Pb') % eksenler isimlendirildi
semilogy(x, Pb_b, 'b', 'LineWidth',1)    % teorik olarak hesaplanan b şıkkı
% çizdirildi
semilogy(x, err_a, 'm:o', 'LineWidth',1) % simulasyon ile elde edilen a 
% şıkkı çizdirildi
semilogy(x, err_b, 'c:o', 'LineWidth',1) % simulasyon ile elde edilen b 
% şıkkı çizdirildi
legend('(a) teorik', '(b) teorik', '(a) simulasyon', '(b) simulasyon')
hold off