function Y_smooth = interpUpsampleSmoothPressure(X, Y, pressure_target, method, smooth_window)
    % INTERPUPSAMPLESMOOTHPRESSURE Interpole, upsample et lisse des courbes de fréquence en fonction de la pression
    %   X : Vecteur des pressions originales (axe x, ex: [0, 1, 2, 3] bars)
    %   Y : Matrice des fréquences (axe y, chaque colonne = un neurone, chaque ligne = une pression)
    %   pressure_target : Vecteur des pressions cibles (ex: 0:0.1:3 pour upsampler)
    %   method : Méthode d'interpolation ('linear', 'spline', 'pchip', etc.)
    %   smooth_window : Taille de la fenêtre de lissage (ex: 3 ou 5)
    %
    %   Y_smooth : Matrice des fréquences interpolées, upsamplées et lissées (taille = length(pressure_target) x size(Y,2))

    % 1. Vérifications
    if isempty(X) || isempty(Y)
        error('X ou Y est vide.');
    end
    if ~isvector(X)
        error('X doit être un vecteur (pression).');
    end
    if size(Y, 1) ~= length(X)
        error('Le nombre de lignes de Y doit correspondre à la longueur de X.');
    end

    % 2. Initialiser la matrice de sortie
    Y_smooth = zeros(length(pressure_target), size(Y, 2));

    % 3. Traiter chaque colonne de Y (chaque neurone)
    for col = 1:size(Y, 2)
        y = Y(:, col);

        % Interpolation
        y_interp = interp1(X, y, pressure_target, method, 'extrap');

        % Lissage (moyenne mobile)
        if exist('smooth_window', 'var') && smooth_window > 1
            window = ones(1, smooth_window) / smooth_window;
            y_smooth_col = conv(y_interp, window, 'same');
        else
            y_smooth_col = y_interp;
        end

        Y_smooth(:, col) = y_smooth_col;
    end

    % 4. Affichage optionnel (pour la première colonne)
    figure;
    subplot(2,1,1);
    plot(X, Y(:, 1), 'o-', 'DisplayName', 'Original');
    hold on;
    plot(pressure_target, interp1(X, Y(:, 1), pressure_target, method), 'x-', 'DisplayName', 'Interpolé');
    title(['Avant lissage - Neurone 1']);
    xlabel('Pression');
    ylabel('Fréquence (Hz)');
    legend;

    subplot(2,1,2);
    plot(pressure_target, Y_smooth(:, 1), 'r-', 'DisplayName', 'Lissé');
    title(['Après lissage - Neurone 1']);
    xlabel('Pression');
    ylabel('Fréquence (Hz)');
    legend;
end