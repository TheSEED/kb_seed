# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.097",
	package_date => 1403125823,
	package_date_str => "Jun 18, 2014 16:10:23",
    };
    return bless $self, $class;
}
1;
