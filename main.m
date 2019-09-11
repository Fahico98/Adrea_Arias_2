%% a. Funcion escalon unitario. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

fprintf('Ingrese los valores del vector de tiempo...\n');

while true
    t_inicial = input('valor inicial: ');
    t_final = input('valor final: ');
    if t_inicial >= t_final
        fprintf('\nEl valor inicial del vector debe ser menor al valor final.\n');
        fprintf('Por faovr intentelo de nuevo...\n');
    else
        break;
    end
end

t = t_inicial : 0.01 : t_final;
u = u(t);

plot(t, u, 'LineWidth', 2, 'Color', 'r');
xlabel("t(s)");
ylabel("y");
grid on;

%% b. Funcion rampa unitaria. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

fprintf('Ingrese los valores del vector de tiempo...\n');

while true
    t_inicial = input('valor inicial: ');
    t_final = input('valor final: ');
    if t_inicial >= t_final
        fprintf('\nEl valor inicial del vector debe ser menor al valor final.\n');
        fprintf('Por faovr intentelo de nuevo...\n');
    else
        break;
    end
end

t = t_inicial : 0.01 : t_final;
ramp = ramp(t);

plot(t, ramp, 'LineWidth', 2, 'Color', 'r');
xlabel("t(s)");
ylabel("y");
grid on;

%% c. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a. g(t) = 0.5 * sen(4*pi*t) * u(t)

clear all;
clc;

dt = 0.01;
t = -5 : dt : 5;
u = u(t);

g = 0.5 .* sin(4*pi*t) .* u;
G = cumtrapz(t, g);
dg = diff(g)/dt;
dg(length(dg) + 1) = 0;

subplot(3, 1, 1);
plot(t, g, 'LineWidth', 2, 'Color', 'k');
xlabel("t(s)");
ylabel("g");
grid on;

subplot(3, 1, 2);
plot(t, dg, 'LineWidth', 2, 'Color', 'r');
xlabel("t(s)");
ylabel("dg");
grid on;

subplot(3, 1, 3);
plot(t, G, 'LineWidth', 2, 'Color', 'b');
xlabel("t(s)");
ylabel("G");
grid on;

E = sum(abs(g) .^ 2);

syms N;
P = limit(E/(2*N + 1), N, inf);

fprintf('Energia de la señal: %g\n', E);
fprintf('Potencia de la señal: %d\n', P);


%% c. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b. g(t) = 2 * u(4 - t)

clear all;
clc;

dt = 0.01;
t = -5 : dt : 5;
u = u(4 - t);

g = 2 .* u;
G = cumtrapz(t, g);
dg = diff(g)/dt;
dg(length(dg) + 1) = 0;

subplot(3, 1, 1);
plot(t, g, 'LineWidth', 2, 'Color', 'k');
xlabel("t(s)");
ylabel("g");
grid on;

subplot(3, 1, 2);
plot(t, dg, 'LineWidth', 2, 'Color', 'r');
xlabel("t(s)");
ylabel("dg");
grid on;

subplot(3, 1, 3);
plot(t, G, 'LineWidth', 2, 'Color', 'b');
xlabel("t(s)");
ylabel("G");
grid on;

E = sum(abs(g) .^ 2);

syms N;
P = limit(E/(2*N + 1), N, inf);

fprintf('Energia de la señal: %g\n', E);
fprintf('Potencia de la señal: %d\n', P);

%% c. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c. g(t) = -4 * ramp(t) * u(t - 2)

clear all;
clc;

dt = 0.01;
t = -5 : dt : 5;
u = u(t - 2);
ramp = ramp(t);

g = -4 .* ramp .* u;
G = cumtrapz(t, g);
dg = diff(g)/dt;
dg(length(dg) + 1) = 0;

subplot(3, 1, 1);
plot(t, g, 'LineWidth', 2, 'Color', 'k');
xlabel("t(s)");
ylabel("g");
grid on;

subplot(3, 1, 2);
plot(t, dg, 'LineWidth', 2, 'Color', 'r');
xlabel("t(s)");
ylabel("dg");
grid on;

subplot(3, 1, 3);
plot(t, G, 'LineWidth', 2, 'Color', 'b');
xlabel("t(s)");
ylabel("G");
grid on;

E = sum(abs(g) .^ 2);

syms N;
P = limit(E/(2*N + 1), N, inf);

fprintf('Energia de la señal: %g\n', E);
fprintf('Potencia de la señal: %d\n', P);

%% c. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d. g(t) = u(t + 0.5) * ramp(0.5 - t)

clear all;
clc;

dt = 0.01;
t = -5 : dt : 5;
u = u(t + 0.5);
ramp = ramp(0.5 - t);

g = u .* ramp;
G = cumtrapz(t, g);
dg = diff(g)/dt;
dg(length(dg) + 1) = 0;

subplot(3, 1, 1);
plot(t, g, 'LineWidth', 2, 'Color', 'k');
xlabel("t(s)");
ylabel("g");
grid on;

subplot(3, 1, 2);
plot(t, dg, 'LineWidth', 2, 'Color', 'r');
xlabel("t(s)");
ylabel("dg");
grid on;

subplot(3, 1, 3);
plot(t, G, 'LineWidth', 2, 'Color', 'b');
xlabel("t(s)");
ylabel("G");
grid on;

E = sum(abs(g) .^ 2);

syms N;
P = limit(E/(2*N + 1), N, inf);

fprintf('Energia de la señal: %g\n', E);
fprintf('Potencia de la señal: %d\n', P);

%% c. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% e. g(t) = 2 * ramp(t/2) * u(t - 2)

clear all;
clc;

dt = 0.01;
t = -5 : dt : 5;
u = u(t - 2);
ramp = ramp(t/2);

g = 2 .* ramp .* u;
G = cumtrapz(t, g);
dg = diff(g)/dt;
dg(length(dg) + 1) = 0;

subplot(3, 1, 1);
plot(t, g, 'LineWidth', 2, 'Color', 'k');
xlabel("t(s)");
ylabel("g");
grid on;

subplot(3, 1, 2);
plot(t, dg, 'LineWidth', 2, 'Color', 'r');
xlabel("t(s)");
ylabel("dg");
grid on;

subplot(3, 1, 3);
plot(t, G, 'LineWidth', 2, 'Color', 'b');
xlabel("t(s)");
ylabel("G");
grid on;

E = sum(abs(g) .^ 2);

syms N;
P = limit(E/(2*N + 1), N, inf);

fprintf('Energia de la señal: %g\n', E);
fprintf('Potencia de la señal: %d\n', P);

%% c. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f. g(t) = 2 * u(t + 1) - u(t - 1) + sen(4*pi*t)

clear all;
clc;

dt = 0.01;
t = -5 : dt : 5;
u_1 = 2 * u(t + 1);
u_2 = u(t - 1);
u = u_1 - u_2;

g = u + sin(4*pi*t);
G = cumtrapz(t, g);
dg = diff(g)/dt;
dg(length(dg) + 1) = 0;

subplot(3, 1, 1);
plot(t, g, 'LineWidth', 2, 'Color', 'k');
xlabel("t(s)");
ylabel("g");
grid on;

subplot(3, 1, 2);
plot(t, dg, 'LineWidth', 2, 'Color', 'r');
xlabel("t(s)");
ylabel("dg");
grid on;

subplot(3, 1, 3);
plot(t, G, 'LineWidth', 2, 'Color', 'b');
xlabel("t(s)");
ylabel("G");
grid on;

E = sum(abs(g) .^ 2);

syms N;
P = limit(E/(2*N + 1), N, inf);

fprintf('Energia de la señal: %g\n', E);
fprintf('Potencia de la señal: %d\n', P);


