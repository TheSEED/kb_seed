# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.141",
	package_date => 1480520612,
	package_date_str => "Nov 30, 2016 09:43:32",
    };
    return bless $self, $class;
}
1;
