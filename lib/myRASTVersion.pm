# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.146",
	package_date => 1480698699,
	package_date_str => "Dec 02, 2016 11:11:39",
    };
    return bless $self, $class;
}
1;
