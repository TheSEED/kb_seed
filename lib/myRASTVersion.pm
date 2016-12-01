# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.145",
	package_date => 1480628255,
	package_date_str => "Dec 01, 2016 15:37:35",
    };
    return bless $self, $class;
}
1;
