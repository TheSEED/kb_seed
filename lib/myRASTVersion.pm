# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.140",
	package_date => 1480450359,
	package_date_str => "Nov 29, 2016 14:12:39",
    };
    return bless $self, $class;
}
1;
