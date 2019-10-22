f = @(x,y) 100*sqrt(abs(y - 0.01*x^2)) + 0.01*abs(x + 10);
num_var = 2;

max_x = -15;
min_x = 15;
max_y = -3;
min_y = 3;

%{
  A populacao deste Algoritmo genetio sera modelada em uma matriz, onde as
cada linha eh um individuo, e as colunas representam, respectivamente, o genotipo, 
o fenotipo, o f(x) e a probabilidade de selecao para cruzamento.

Cada individuo possui um genotipo binario representando suas variaveis:
10101 10100 10010, onde cada 5 bits sao uma variavel.

O fenotipo eh o equivelente decimal deste binario.

f(x) eh o resultado da aplicacao das variaveis do genotipo na funcao otimizada.

A probabilidade de cruzamento define a chance que um individuo tem de ser escolhido
para cruzar com outros.

Alem disso, sera separado uma Elitizacao de um individuo, ou seja, o 
melhore individuo permanecera inalterado para a proxima geracao.
%}

num_individuos = 10;

fenotipos = repmat(' ', 1, num_individuos);
genotipos = [];
funcs = [];
probs = [];
elite = 1;

% Plota a funcao
X = linspace(min_x, max_x);
Y = linspace(min_y, max_y);

fig = figure();

[xplot,yplot] = meshgrid(X,Y);

zplot = 100*sqrt(abs(yplot - 0.01*xplot.^2)) + 0.01*abs(xplot + 10);

surf(xplot, yplot, zplot);

shading interp;

hold on;
% Cria a primeira geracao de individuos
soma = 0;
for i=1:num_individuos  
  % Cria um fenotipo aleatorio, como uma string
  for b=1: num_var*5
    gene = num2str(randi(2)-1);
    fenotipos(i,b) = gene;
  end
  
  % Cria o genotipo de cada individuo
  genotipos(i) = bin2dec(fenotipos(i,:));

  % Separa as variaveis e calcula a funcao
  variaveis = [];
  for k=1:num_var  
    var = repmat(" ", 1,5);
    % Separa o fenotipo em pedacos de 5 bits
    for j=1 + 5*(k-1):5*k
      var(j-5*(k-1)) = fenotipos(i,j);    
    end
    var = bin2dec(var);
    variaveis(k) = var;
  end
  
  % Calcula f(x)
  fx = f(variaveis(1), variaveis(2));
  funcs(i) = fx;
  soma += fx;
  
  % Encontra a elite de cada geracao
  if fx < funcs(elite)
    elite = i;
  end
end

contador = 0;
while contador < 10
  aux = 1;
  roleta = [];
  for i=1 : num_individuos
    % Calcula a probabilidade de reproducao de cada individuo
    probs(i) = 1 - (soma - funcs(i))/soma;
  
    % Cria a roleta de probabilidades
    for a=aux:round(aux+probs(i)*100) - 1
      roleta(a) = i;
      aux = aux + 1;
    end
  end
  
  break;
  contador = contador + 1;
end

disp(roleta);

