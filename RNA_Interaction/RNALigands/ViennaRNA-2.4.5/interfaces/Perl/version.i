%{
/** @file version.i.in
 * @brief Set $RNA::VERSION to the bindings version
 */
%}

%perlcode %{
our $VERSION = '2.4.5';
sub VERSION () { $VERSION };
%}

