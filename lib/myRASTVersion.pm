# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.143",
	package_date => 1480530109,
	package_date_str => "Nov 30, 2016 12:21:49",
    };
    return bless $self, $class;
}
1;
