function y = getWrongEntries( felix, mirko )
size_felix = size( felix );
size_mirko = size( mirko );
numberOfWrongEntries = 0;
if ( size_felix == size_mirko )
    for i = 1:size_felix(1)
        for j = 1:size_felix(2)
            tmp = felix( i, j ) - mirko( i, j );
            if ( abs( tmp > 1e-10 ) )
                numberOfWrongEntries = numberOfWrongEntries + 1;
                y( numberOfWrongEntries, 1 ) = i;
                y( numberOfWrongEntries, 2 ) = j;
                fprintf( 1, '%2d, %2d: %6f\n', i - 1, j - 1, tmp );
            end
        end
    end
end
y = 0;
