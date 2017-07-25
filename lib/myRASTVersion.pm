# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.151",
	package_date => 1500996992,
	package_date_str => "Jul 25, 2017 10:36:32",
    };
    return bless $self, $class;
}
1;
