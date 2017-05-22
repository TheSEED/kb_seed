# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.150",
	package_date => 1495470100,
	package_date_str => "May 22, 2017 11:21:40",
    };
    return bless $self, $class;
}
1;
