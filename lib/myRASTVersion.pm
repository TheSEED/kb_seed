# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.147",
	package_date => 1480705067,
	package_date_str => "Dec 02, 2016 12:57:47",
    };
    return bless $self, $class;
}
1;
